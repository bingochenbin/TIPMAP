"""Extract REF and matched true ALT sequences for PanPop INS/DEL SV records."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import os
from pathlib import Path
from typing import Iterable, Iterator, List

from tipmap.lib.fasta import (
    IndexedFasta,
    extracted_sequence_to_fasta_record,
    read_fasta,
    write_fasta,
    write_one_fasta_record,
)
from tipmap.lib.matcher import (
    MatchCriteria,
    PairwiseRecordIndex,
    extract_alt_sequence,
    extract_ref_sequence,
    find_matching_alt_records,
    is_supported_tip_sv,
    iter_extracted_sequences,
    should_extract_panpop_alt,
    should_extract_panpop_ref,
)
from tipmap.lib.models import ExtractedSequence, SVRecord
from tipmap.lib.parser import parse_pairwise_vcf, parse_panpop_vcf


@dataclass(frozen=True)
class PairwiseVcfInput:
    """One pairwise genome-to-reference VCF input."""

    genome: str
    path: Path
    query_fasta: Path | None = None


_WORKER_PAIRWISE_RECORDS: List[SVRecord] = []
_WORKER_PAIRWISE_INDEX = PairwiseRecordIndex([])
_WORKER_QUERY_SEQUENCES_BY_GENOME: dict[str, dict[str, str] | IndexedFasta] = {}
_WORKER_REFERENCE_SEQUENCES: dict[str, str] | IndexedFasta = {}
_WORKER_CRITERIA = MatchCriteria()
_WORKER_FLANK = 0
_WORKER_ALT_FLANK = 0
_WORKER_MIN_PANPOP_SEQUENCE_LENGTH = 0
_WORKER_REFERENCE_GENOME = "reference"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panpop-vcf", required=True, help="PanPop SV VCF.")
    parser.add_argument(
        "--pairwise-vcf",
        action="append",
        default=[],
        help="Pairwise genome-to-reference VCF. May be provided multiple times.",
    )
    parser.add_argument(
        "--pairwise-vcf-list",
        help=(
            "Two-column text file listing pairwise VCF inputs as "
            "'genome<TAB>vcf_path', or a three-column file with "
            "'genome<TAB>vcf_path<TAB>query_fasta_path'. Blank lines and lines "
            "beginning with # are ignored. Relative VCF and FASTA paths are "
            "resolved relative to the list file."
        ),
    )
    parser.add_argument("--reference-fasta", required=True, help="Reference genome FASTA.")
    parser.add_argument("--output", required=True, help="Output FASTA path.")
    parser.add_argument("--genome", help="PanPop genome/sample name override.")
    parser.add_argument("--reference-genome", default="reference", help="Reference genome name.")
    parser.add_argument(
        "--max-distance",
        type=int,
        default=100,
        help="Maximum POS/END distance for matching PanPop and pairwise SV records.",
    )
    parser.add_argument(
        "--max-length-ratio-difference",
        type=float,
        default=0.25,
        help="Maximum relative SVLEN difference for matching records.",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=50,
        help="Reference-side flank size to include around the PanPop interval.",
    )
    parser.add_argument(
        "--alt-flank",
        type=int,
        default=50,
        help=(
            "Query-side flank size to include around matched pairwise ALT intervals. "
            "Requires ASM_* fields in pairwise VCF records and query_fasta_path in "
            "--pairwise-vcf-list. Records without these inputs fall back to VCF ALT."
        ),
    )
    parser.add_argument(
        "--min-panpop-sequence-length",
        type=int,
        default=50,
        help=(
            "Minimum original PanPop REF/ALT allele length required for sequence "
            "output. Filtering is role-specific and is applied before flank "
            "extension: short REF alleles do not block long ALT extraction, and "
            "short ALT alleles do not block long REF extraction."
        ),
    )
    parser.add_argument(
        "--fasta-index-mode",
        choices=("auto", "require", "rebuild"),
        default="auto",
        help=(
            "FASTA index handling for reference and Query FASTA files. 'auto' uses "
            "or builds .fai indexes, 'require' fails if an index is missing, and "
            "'rebuild' regenerates indexes."
        ),
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=100,
        help="Number of PanPop records per multiprocessing chunk.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help=(
            "Number of worker processes for pairwise VCF parsing and per-SV "
            "sequence extraction. Use 1 for serial execution."
        ),
    )
    return parser


def read_pairwise_vcf_list(path: str | Path) -> List[PairwiseVcfInput]:
    """Read a genome/VCF/query-FASTA manifest for pairwise VCF inputs."""

    manifest_path = Path(path)
    inputs: List[PairwiseVcfInput] = []
    with manifest_path.open("rt", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = _split_manifest_line(line)
            if fields is None:
                continue
            genome, vcf_path_raw, query_fasta_raw = fields
            if not genome:
                msg = "Missing genome name in %s line %d" % (manifest_path, line_number)
                raise ValueError(msg)
            if not vcf_path_raw:
                msg = "Missing VCF path in %s line %d" % (manifest_path, line_number)
                raise ValueError(msg)
            vcf_path = Path(vcf_path_raw)
            if not vcf_path.is_absolute():
                vcf_path = manifest_path.parent / vcf_path
            query_fasta = None
            if query_fasta_raw:
                query_fasta = Path(query_fasta_raw)
                if not query_fasta.is_absolute():
                    query_fasta = manifest_path.parent / query_fasta
            inputs.append(PairwiseVcfInput(genome=genome, path=vcf_path, query_fasta=query_fasta))
    if not inputs:
        msg = "Pairwise VCF list contains no inputs: %s" % manifest_path
        raise ValueError(msg)
    return inputs


def collect_pairwise_records(
    direct_vcf_paths: Iterable[str],
    *,
    pairwise_vcf_list: str | None = None,
    workers: int = 1,
) -> List[SVRecord]:
    """Parse pairwise VCF records from direct paths and an optional manifest."""

    inputs = _collect_pairwise_inputs(direct_vcf_paths, pairwise_vcf_list=pairwise_vcf_list)
    return _collect_pairwise_records_from_inputs(inputs, workers=workers)


def _collect_pairwise_records_from_inputs(
    inputs: List[PairwiseVcfInput],
    *,
    workers: int = 1,
) -> List[SVRecord]:
    """Parse pairwise VCF records from collected input descriptors."""

    if not inputs:
        return []
    if workers < 1:
        msg = "workers must be >= 1"
        raise ValueError(msg)
    if workers == 1 or len(inputs) == 1:
        nested_records = [_parse_pairwise_input(input_record) for input_record in inputs]
    else:
        max_workers = min(workers, len(inputs), os.cpu_count() or workers)
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            nested_records = list(executor.map(_parse_pairwise_input, inputs))

    pairwise_records = []
    for records in nested_records:
        pairwise_records.extend(records)
    return pairwise_records


def _collect_pairwise_inputs(
    direct_vcf_paths: Iterable[str],
    *,
    pairwise_vcf_list: str | None = None,
) -> List[PairwiseVcfInput]:
    inputs: List[PairwiseVcfInput] = []
    for path in direct_vcf_paths:
        inputs.append(PairwiseVcfInput(genome="", path=Path(path)))
    if pairwise_vcf_list is not None:
        inputs.extend(read_pairwise_vcf_list(pairwise_vcf_list))
    return inputs


def _parse_pairwise_input(input_record: PairwiseVcfInput) -> List[SVRecord]:
    genome = input_record.genome or None
    return list(parse_pairwise_vcf(input_record.path, genome=genome))


def collect_query_sequences_by_genome(
    pairwise_inputs: Iterable[PairwiseVcfInput],
) -> dict[str, dict[str, str]]:
    """Read Query genome FASTA files listed beside pairwise VCF inputs."""

    query_sequences_by_genome: dict[str, dict[str, str]] = {}
    for input_record in pairwise_inputs:
        if input_record.query_fasta is None or not input_record.genome:
            continue
        query_sequences_by_genome[input_record.genome] = read_fasta(input_record.query_fasta)
    return query_sequences_by_genome


def collect_query_indexes_by_genome(
    pairwise_inputs: Iterable[PairwiseVcfInput],
    *,
    index_mode: str = "auto",
) -> dict[str, IndexedFasta]:
    """Open Query genome FASTA indexes listed beside pairwise VCF inputs."""

    query_indexes_by_genome: dict[str, IndexedFasta] = {}
    for input_record in pairwise_inputs:
        if input_record.query_fasta is None or not input_record.genome:
            continue
        query_indexes_by_genome[input_record.genome] = IndexedFasta.open(
            input_record.query_fasta,
            index_mode=index_mode,
        )
    return query_indexes_by_genome


def iter_collected_extracted_sequences(
    panpop_records: Iterable[SVRecord],
    pairwise_records: Iterable[SVRecord],
    reference_sequences: dict[str, str] | IndexedFasta,
    *,
    criteria: MatchCriteria,
    flank: int = 0,
    alt_flank: int = 0,
    min_panpop_sequence_length: int = 0,
    query_sequences_by_genome: dict[str, dict[str, str] | IndexedFasta] | None = None,
    reference_genome: str = "reference",
    workers: int = 1,
    chunk_size: int = 100,
) -> Iterator[ExtractedSequence]:
    """Yield REF/ALT sequences serially or by PanPop SV in worker processes."""

    panpop_records_list = list(panpop_records)
    pairwise_records_list = list(pairwise_records)
    if workers < 1:
        msg = "workers must be >= 1"
        raise ValueError(msg)
    if chunk_size < 1:
        msg = "chunk_size must be >= 1"
        raise ValueError(msg)
    if alt_flank < 0:
        msg = "alt_flank must be >= 0"
        raise ValueError(msg)
    if min_panpop_sequence_length < 0:
        msg = "min_panpop_sequence_length must be >= 0"
        raise ValueError(msg)

    query_sequences = query_sequences_by_genome or {}
    if workers == 1 or len(panpop_records_list) <= 1:
        pairwise_index = PairwiseRecordIndex(pairwise_records_list)
        yield from iter_extracted_sequences(
            panpop_records_list,
            pairwise_index,
            reference_sequences,
            criteria=criteria,
            flank=flank,
            alt_flank=alt_flank,
            min_panpop_sequence_length=min_panpop_sequence_length,
            query_sequences_by_genome=query_sequences,
            reference_genome=reference_genome,
        )
        return

    max_workers = min(workers, len(panpop_records_list), os.cpu_count() or workers)
    with ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=_init_extraction_worker,
        initargs=(
            pairwise_records_list,
            reference_sequences,
            query_sequences,
            criteria,
            flank,
            alt_flank,
            min_panpop_sequence_length,
            reference_genome,
        ),
    ) as executor:
        nested_records = executor.map(
            _extract_sequences_for_panpop,
            panpop_records_list,
            chunksize=chunk_size,
        )
        for records in nested_records:
            yield from records


def collect_extracted_sequences(
    panpop_records: Iterable[SVRecord],
    pairwise_records: Iterable[SVRecord],
    reference_sequences: dict[str, str] | IndexedFasta,
    *,
    criteria: MatchCriteria,
    flank: int = 0,
    alt_flank: int = 0,
    min_panpop_sequence_length: int = 0,
    query_sequences_by_genome: dict[str, dict[str, str] | IndexedFasta] | None = None,
    reference_genome: str = "reference",
    workers: int = 1,
    chunk_size: int = 100,
) -> List[ExtractedSequence]:
    """Collect extracted sequences into a list for tests and small runs."""

    return list(
        iter_collected_extracted_sequences(
            panpop_records,
            pairwise_records,
            reference_sequences,
            criteria=criteria,
            flank=flank,
            alt_flank=alt_flank,
            min_panpop_sequence_length=min_panpop_sequence_length,
            query_sequences_by_genome=query_sequences_by_genome,
            reference_genome=reference_genome,
            workers=workers,
            chunk_size=chunk_size,
        )
    )


def _init_extraction_worker(
    pairwise_records: List[SVRecord],
    reference_sequences: dict[str, str] | IndexedFasta,
    query_sequences_by_genome: dict[str, dict[str, str] | IndexedFasta],
    criteria: MatchCriteria,
    flank: int,
    alt_flank: int,
    min_panpop_sequence_length: int,
    reference_genome: str,
) -> None:
    global _WORKER_PAIRWISE_RECORDS
    global _WORKER_PAIRWISE_INDEX
    global _WORKER_QUERY_SEQUENCES_BY_GENOME
    global _WORKER_REFERENCE_SEQUENCES
    global _WORKER_CRITERIA
    global _WORKER_FLANK
    global _WORKER_ALT_FLANK
    global _WORKER_MIN_PANPOP_SEQUENCE_LENGTH
    global _WORKER_REFERENCE_GENOME

    _WORKER_PAIRWISE_RECORDS = pairwise_records
    _WORKER_PAIRWISE_INDEX = PairwiseRecordIndex(pairwise_records)
    _WORKER_QUERY_SEQUENCES_BY_GENOME = query_sequences_by_genome
    _WORKER_REFERENCE_SEQUENCES = reference_sequences
    _WORKER_CRITERIA = criteria
    _WORKER_FLANK = flank
    _WORKER_ALT_FLANK = alt_flank
    _WORKER_MIN_PANPOP_SEQUENCE_LENGTH = min_panpop_sequence_length
    _WORKER_REFERENCE_GENOME = reference_genome


def _extract_sequences_for_panpop(panpop_record: SVRecord) -> List[ExtractedSequence]:
    if not is_supported_tip_sv(panpop_record):
        return []

    records: List[ExtractedSequence] = []
    if should_extract_panpop_ref(panpop_record, _WORKER_MIN_PANPOP_SEQUENCE_LENGTH):
        records.append(
            extract_ref_sequence(
                panpop_record,
                _WORKER_REFERENCE_SEQUENCES,
                flank=_WORKER_FLANK,
                reference_genome=_WORKER_REFERENCE_GENOME,
            )
        )
    if not should_extract_panpop_alt(panpop_record, _WORKER_MIN_PANPOP_SEQUENCE_LENGTH):
        return records
    for pairwise_record in find_matching_alt_records(
        panpop_record,
        _WORKER_PAIRWISE_INDEX,
        criteria=_WORKER_CRITERIA,
    ):
        records.append(
            extract_alt_sequence(
                panpop_record,
                pairwise_record,
                query_sequences_by_genome=_WORKER_QUERY_SEQUENCES_BY_GENOME,
                alt_flank=_WORKER_ALT_FLANK,
            )
        )
    return records


def _split_manifest_line(line: str) -> tuple[str, str, str | None] | None:
    if "\t" in line:
        fields = [field.strip() for field in line.split("\t")]
    else:
        fields = line.split()
    if len(fields) not in {2, 3}:
        msg = "Pairwise VCF list lines must contain two or three columns"
        raise ValueError(msg)
    genome, path = fields[:2]
    query_fasta = fields[2] if len(fields) == 3 else None
    if genome.lower() in {"genome", "sample", "sample_id"} and path.lower() in {
        "vcf",
        "path",
        "vcf_path",
    }:
        return None
    return genome, path, query_fasta


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()
    if not args.pairwise_vcf and args.pairwise_vcf_list is None:
        msg = "At least one --pairwise-vcf or --pairwise-vcf-list is required."
        parser.error(msg)

    reference_sequences = IndexedFasta.open(args.reference_fasta, index_mode=args.fasta_index_mode)
    panpop_records = list(parse_panpop_vcf(args.panpop_vcf, genome=args.genome))
    pairwise_inputs = _collect_pairwise_inputs(
        args.pairwise_vcf,
        pairwise_vcf_list=args.pairwise_vcf_list,
    )
    pairwise_records = _collect_pairwise_records_from_inputs(pairwise_inputs, workers=args.workers)
    query_sequences_by_genome: dict[str, IndexedFasta] = {}
    if args.alt_flank > 0:
        query_sequences_by_genome = collect_query_indexes_by_genome(
            pairwise_inputs,
            index_mode=args.fasta_index_mode,
        )

    criteria = MatchCriteria(
        max_distance=args.max_distance,
        max_length_ratio_difference=args.max_length_ratio_difference,
    )
    ref_count = 0
    alt_count = 0
    with Path(args.output).open("wt", encoding="utf-8", newline="\n") as handle:
        for record in iter_collected_extracted_sequences(
            panpop_records,
            pairwise_records,
            reference_sequences,
            criteria=criteria,
            flank=args.flank,
            alt_flank=args.alt_flank,
            min_panpop_sequence_length=args.min_panpop_sequence_length,
            query_sequences_by_genome=query_sequences_by_genome,
            reference_genome=args.reference_genome,
            workers=args.workers,
            chunk_size=args.chunk_size,
        ):
            if record.role == "ref":
                ref_count += 1
            elif record.role == "alt":
                alt_count += 1
            write_one_fasta_record(extracted_sequence_to_fasta_record(record), handle)

    print(
        "Extracted %d REF sequences and %d true ALT sequences to %s."
        % (ref_count, alt_count, args.output)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
