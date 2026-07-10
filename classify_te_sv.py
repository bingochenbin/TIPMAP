"""Classify PanPop SV alleles as TE-derived and write a TIP-map VCF."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path
import re
import shlex
import subprocess
from typing import Callable, Iterable, Iterator, Sequence, TextIO

from tipmap.lib.fasta import FastaRecord, iter_fasta, read_fasta, reverse_complement, write_fasta_records
from tipmap.lib.matcher import panpop_sv_key
from tipmap.lib.parser import parse_panpop_vcf, parse_vcf_record_line
from tipmap.lib.utils import sequence_md5

TIP_INFO_LINES = [
    '##INFO=<ID=TIPMAP,Number=0,Type=Flag,Description="Record retained by TIPMap as a TE-SV/TIP candidate">',
    '##INFO=<ID=TIP_TE_REF,Number=1,Type=Integer,Description="1 if the REF allele is classified as TE-derived">',
    '##INFO=<ID=TIP_TE_ALTS,Number=.,Type=Integer,Description="Original 1-based ALT allele numbers classified as TE-derived">',
    '##INFO=<ID=TIP_RETAINED_ALTS,Number=.,Type=Integer,Description="Original 1-based ALT allele numbers retained in this TIP VCF record">',
]
BLAST_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
BLAST_OUTFMT6_COLUMNS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)
CommandRunner = Callable[[Sequence[str], Path | None], subprocess.CompletedProcess[str]]


@dataclass(frozen=True)
class TEMetadata:
    seq_id: str = ""
    family: str = ""
    te_class: str = ""
    start: str = ""
    end: str = ""
    strand: str = ""
    feature_type: str = ""
    attributes: str = ""


@dataclass(frozen=True)
class TEHit:
    start: int
    end: int
    identity: float
    family: str = ""
    te_class: str = ""
    subject_id: str = ""
    metadata: TEMetadata | None = None


@dataclass(frozen=True)
class AlleleEvidence:
    md5: str
    allele_length: int
    te_covered_bp: int
    coverage: float
    weighted_identity: float | None
    is_te: bool
    family: str = ""
    te_class: str = ""
    supporting_te_annotations: str = ""


@dataclass(frozen=True)
class AlleleReportRow:
    record_id: str
    chrom: str
    pos: str
    svtype: str
    allele_role: str
    original_allele_number: int
    new_allele_number: str
    allele_length: int
    md5: str
    te_covered_bp: int
    te_coverage: float
    weighted_identity: str
    is_te: bool
    retained: bool
    family: str
    te_class: str
    supporting_te_annotations: str


@dataclass(frozen=True)
class ProcessedVcfRecord:
    line: str | None
    reports: list[AlleleReportRow] = field(default_factory=list)


@dataclass(frozen=True)
class PanpopAlleleIndexEntry:
    md5: str
    record_id: str
    chrom: str
    pos: int
    end: int
    svtype: str
    allele_role: str
    allele_number: int
    length: int


@dataclass
class DedupTEFragmentEntry:
    md5: str
    representative_seq_id: str
    family: str
    te_class: str
    all_families: set[str] = field(default_factory=set)
    all_classes: set[str] = field(default_factory=set)
    source_count: int = 0


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panpop-vcf", required=True, help="Input PanPop merged SV VCF.")
    parser.add_argument(
        "--sv-fasta",
        help="SV FASTA annotated by EDTA/TE tools. Required unless --blast-tsv is provided.",
    )
    parser.add_argument(
        "--te-annotations",
        help=(
            "TE annotation TSV from annotate_te.py or scripts/parse_edta_gff3.py. "
            "Required for integrated BLAST; optional with --blast-tsv for family metadata."
        ),
    )
    parser.add_argument(
        "--blast-tsv",
        help="Existing BLASTN outfmt 6 result. If omitted, classify_te_sv.py runs BLASTN internally.",
    )
    parser.add_argument("--output-vcf", required=True, help="Output TIP-map VCF path.")
    parser.add_argument("--allele-report", help="Optional allele-level TE evidence TSV path.")
    parser.add_argument("--workdir", help="Working directory for allele FASTA, TE FASTA, and BLAST output.")
    parser.add_argument("--panpop-alleles-fasta", help="Optional path for generated PanPop allele FASTA.")
    parser.add_argument("--te-fragments-fasta", help="Optional path for generated TE fragment FASTA.")
    parser.add_argument("--makeblastdb", default="makeblastdb", help="makeblastdb executable path. Default: makeblastdb.")
    parser.add_argument("--blast-db-prefix", help="Optional BLAST database prefix. Default: <workdir>/te_fragments_db.")
    parser.add_argument("--blastn", default="blastn", help="blastn executable path. Default: blastn.")
    parser.add_argument("--blast-threads", type=int, default=1, help="Threads passed to blastn -num_threads. Default: 1.")
    parser.add_argument("--blast-evalue", type=float, default=1e-5, help="Maximum BLASTN e-value. Default: 1e-5.")
    parser.add_argument(
        "--blast-arg",
        action="append",
        default=[],
        help="Extra argument string passed to blastn. Repeat for multiple argument groups.",
    )
    parser.add_argument("--min-panpop-allele-length", type=int, default=50, help="Minimum PanPop allele length for BLAST. Default: 50.")
    parser.add_argument("--min-te-fragment-length", type=int, default=10, help="Minimum extracted TE fragment length. Default: 10.")
    parser.add_argument("--fasta-line-width", type=int, default=80, help="Generated FASTA line width. Default: 80.")
    parser.add_argument(
        "--min-te-coverage",
        type=float,
        default=0.60,
        help="Minimum union TE coverage over the PanPop allele sequence. Default: 0.60.",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=80.0,
        help="Minimum BLASTN weighted identity in percent. Default: 80.",
    )
    parser.add_argument(
        "--min-te-covered-bp",
        type=int,
        default=40,
        help="Minimum union TE-covered base pairs. Default: 40.",
    )
    return parser


def iter_panpop_te_alleles(panpop_vcf: str | Path, *, min_length: int = 1) -> Iterator[FastaRecord]:
    """Yield PanPop alleles for TE classification: INS ALT, DEL REF."""

    for record, _entry in iter_panpop_te_alleles_with_index(panpop_vcf, min_length=min_length):
        yield record


def iter_panpop_te_alleles_with_index(
    panpop_vcf: str | Path,
    *,
    min_length: int = 1,
) -> Iterator[tuple[FastaRecord, PanpopAlleleIndexEntry]]:
    """Yield PanPop allele FASTA records with reusable MD5 index entries."""

    if min_length < 1:
        msg = "min_length must be >= 1"
        raise ValueError(msg)
    emitted_del_refs: set[tuple[str, int, str, str]] = set()
    for record in parse_panpop_vcf(panpop_vcf):
        sv_id = panpop_sv_key(record)
        if record.svtype == "INS":
            sequence = record.alt
            role = "alt"
            original_allele = record.allele_index + 1
        elif record.svtype == "DEL":
            sequence = record.ref
            role = "ref"
            original_allele = 0
            del_key = (record.chrom, record.pos, record.source_id or sv_id, sequence_md5(sequence))
            if del_key in emitted_del_refs:
                continue
            emitted_del_refs.add(del_key)
        else:
            continue
        if len(sequence) < min_length:
            continue
        md5 = sequence_md5(sequence)
        name = "%s|%s|panpop|%s|%d|%d|%s|%s" % (
            sv_id,
            role,
            record.chrom,
            record.pos,
            record.end,
            record.svtype,
            md5,
        )
        description = "source_id=%s allele_index=%d original_allele=%d length=%d" % (
            record.source_id or ".",
            record.allele_index,
            original_allele,
            len(sequence),
        )
        index_entry = PanpopAlleleIndexEntry(
            md5=md5,
            record_id=record.source_id or sv_id,
            chrom=record.chrom,
            pos=record.pos,
            end=record.end,
            svtype=record.svtype,
            allele_role=role,
            allele_number=original_allele,
            length=len(sequence),
        )
        yield FastaRecord(name=name, sequence=sequence, description=description), index_entry


def write_panpop_te_alleles(
    panpop_vcf: str | Path,
    output: str | Path,
    *,
    min_length: int = 1,
    line_width: int = 80,
) -> int:
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with output_path.open("wt", encoding="utf-8", newline="\n") as handle:
        for record in iter_panpop_te_alleles(panpop_vcf, min_length=min_length):
            write_fasta_records([record], handle, line_width=line_width)
            count += 1
    return count


def write_panpop_te_alleles_and_index(
    panpop_vcf: str | Path,
    output: str | Path,
    index_output: str | Path,
    *,
    min_length: int = 1,
    line_width: int = 80,
) -> dict[tuple[str, int, int, str, str, int], PanpopAlleleIndexEntry]:
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    entries: list[PanpopAlleleIndexEntry] = []
    lookup: dict[tuple[str, int, int, str, str, int], PanpopAlleleIndexEntry] = {}
    with output_path.open("wt", encoding="utf-8", newline="\n") as handle:
        for record, entry in iter_panpop_te_alleles_with_index(panpop_vcf, min_length=min_length):
            write_fasta_records([record], handle, line_width=line_width)
            entries.append(entry)
            lookup[panpop_allele_index_key(entry.chrom, entry.pos, entry.end, entry.svtype, entry.allele_role, entry.allele_number)] = entry
    write_panpop_allele_index(entries, index_output)
    return lookup


def write_panpop_allele_index(entries: Sequence[PanpopAlleleIndexEntry], output: str | Path) -> None:
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["md5", "record_id", "chrom", "pos", "end", "svtype", "allele_role", "allele_number", "length"],
        )
        writer.writeheader()
        for entry in entries:
            writer.writerow(entry.__dict__)


def panpop_allele_index_key(
    chrom: str,
    pos: int,
    end: int,
    svtype: str,
    allele_role: str,
    allele_number: int,
) -> tuple[str, int, int, str, str, int]:
    return (chrom, pos, end, svtype, allele_role, allele_number)


def iter_te_fragments(
    sv_fasta: str | Path,
    te_annotations: str | Path,
    *,
    min_length: int = 1,
) -> Iterator[FastaRecord]:
    """Yield TE fragments from an annotated SV FASTA and TE annotation TSV."""

    if min_length < 1:
        msg = "min_length must be >= 1"
        raise ValueError(msg)
    sequences = read_fasta(sv_fasta)
    with Path(te_annotations).open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return
        for row in reader:
            seq_id = row.get("seq_id", "")
            if not seq_id or seq_id not in sequences:
                continue
            start = _parse_int(row.get("start"))
            end = _parse_int(row.get("end"))
            if start is None or end is None:
                continue
            sequence = sequences[seq_id]
            left = max(min(start, end), 1)
            right = min(max(start, end), len(sequence))
            if right < left:
                continue
            fragment = sequence[left - 1 : right]
            if row.get("strand") == "-":
                fragment = reverse_complement(fragment)
            if len(fragment) < min_length:
                continue
            fragment_md5 = sequence_md5(fragment)
            name = "%s::te:%d-%d:%s" % (seq_id, left, right, fragment_md5)
            description = "family=%s class=%s strand=%s length=%d" % (
                row.get("family", "") or ".",
                row.get("class", "") or ".",
                row.get("strand", ".") or ".",
                len(fragment),
            )
            yield FastaRecord(name=name, sequence=fragment, description=description)


def write_te_fragments(
    sv_fasta: str | Path,
    te_annotations: str | Path,
    output: str | Path,
    *,
    min_length: int = 1,
    line_width: int = 80,
) -> int:
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with output_path.open("wt", encoding="utf-8", newline="\n") as handle:
        for record in iter_te_fragments(sv_fasta, te_annotations, min_length=min_length):
            write_fasta_records([record], handle, line_width=line_width)
            count += 1
    return count


def write_deduplicated_te_fragments(
    input_fasta: str | Path,
    output_fasta: str | Path,
    metadata_output: str | Path,
    *,
    line_width: int = 80,
) -> int:
    """Write sequence-level de-duplicated TE fragments and classification metadata."""

    output_path = Path(output_fasta)
    metadata_path = Path(metadata_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    records_by_md5: dict[str, FastaRecord] = {}
    metadata_by_md5: dict[str, DedupTEFragmentEntry] = {}
    for record in iter_fasta(input_fasta):
        md5 = sequence_md5(record.sequence)
        family, te_class = parse_te_fragment_description(record.description or "")
        if md5 not in records_by_md5:
            records_by_md5[md5] = record
            metadata_by_md5[md5] = DedupTEFragmentEntry(
                md5=md5,
                representative_seq_id=record.name,
                family=family,
                te_class=te_class,
                all_families=set(),
                all_classes=set(),
                source_count=0,
            )
        entry = metadata_by_md5[md5]
        if family:
            entry.all_families.add(family)
        if te_class:
            entry.all_classes.add(te_class)
        entry.source_count += 1
    with output_path.open("wt", encoding="utf-8", newline="\n") as handle:
        for md5 in sorted(records_by_md5):
            write_fasta_records([records_by_md5[md5]], handle, line_width=line_width)
    write_deduplicated_te_fragment_metadata(metadata_by_md5.values(), metadata_path)
    return len(records_by_md5)


def parse_te_fragment_description(description: str) -> tuple[str, str]:
    values: dict[str, str] = {}
    for item in description.split():
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        values[key] = "" if value == "." else value
    return values.get("family", ""), values.get("class", "")


def summarize_te_label(values: set[str], first_value: str) -> str:
    cleaned = sorted(value for value in values if value)
    if not cleaned:
        return ""
    if len(cleaned) == 1:
        return cleaned[0]
    return "mixed" if first_value not in cleaned or len(cleaned) > 1 else first_value


def write_deduplicated_te_fragment_metadata(entries: Iterable[DedupTEFragmentEntry], output: str | Path) -> None:
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["md5", "representative_seq_id", "family", "class", "all_families", "all_classes", "source_count"],
        )
        writer.writeheader()
        for entry in sorted(entries, key=lambda item: item.md5):
            all_families = sorted(value for value in entry.all_families if value)
            all_classes = sorted(value for value in entry.all_classes if value)
            writer.writerow(
                {
                    "md5": entry.md5,
                    "representative_seq_id": entry.representative_seq_id,
                    "family": summarize_te_label(entry.all_families, entry.family),
                    "class": summarize_te_label(entry.all_classes, entry.te_class),
                    "all_families": ",".join(all_families),
                    "all_classes": ",".join(all_classes),
                    "source_count": entry.source_count,
                }
            )


def run_makeblastdb(
    *,
    fasta: str | Path,
    db_prefix: str | Path,
    makeblastdb: str = "makeblastdb",
    runner: CommandRunner | None = None,
) -> Path:
    db_path = Path(db_prefix)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    command = [makeblastdb, "-in", str(fasta), "-dbtype", "nucl", "-out", str(db_path)]
    result = (runner or _default_runner)(command, None)
    if result.returncode != 0:
        msg = "makeblastdb failed with exit code %d" % result.returncode
        raise RuntimeError(msg)
    return db_path


def run_blastn(
    *,
    query_fasta: str | Path,
    blast_db: str | Path,
    output: str | Path,
    blastn: str = "blastn",
    blast_threads: int = 1,
    blast_evalue: float = 1e-5,
    extra_args: Sequence[str] = (),
    runner: CommandRunner | None = None,
) -> Path:
    if blast_threads < 1:
        msg = "blast_threads must be >= 1"
        raise ValueError(msg)
    if blast_evalue <= 0:
        msg = "blast_evalue must be > 0"
        raise ValueError(msg)
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        blastn,
        "-query",
        str(query_fasta),
        "-db",
        str(blast_db),
        "-out",
        str(output_path),
        "-outfmt",
        BLAST_OUTFMT,
        "-evalue",
        str(blast_evalue),
        "-num_threads",
        str(blast_threads),
    ]
    for extra in extra_args:
        command.extend(shlex.split(extra))
    result = (runner or _default_runner)(command, None)
    if result.returncode != 0:
        msg = "blastn failed with exit code %d" % result.returncode
        raise RuntimeError(msg)
    return output_path


def read_te_metadata(path: str | Path | None) -> dict[str, TEMetadata]:
    """Read optional TE annotation metadata keyed by sequence ID and MD5."""

    if path is None:
        return {}
    metadata: dict[str, TEMetadata] = {}
    with Path(path).open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return metadata
        for row in reader:
            seq_id = row.get("seq_id", "")
            md5 = row.get("md5", "") or md5_from_identifier(seq_id)
            item = TEMetadata(
                seq_id=seq_id,
                family=row.get("family", ""),
                te_class=row.get("class", ""),
                start=row.get("start", ""),
                end=row.get("end", ""),
                strand=row.get("strand", ""),
                feature_type=row.get("type", "") or row.get("feature_type", ""),
                attributes=row.get("attributes", ""),
            )
            if seq_id:
                metadata[seq_id] = item
                fragment_key = te_fragment_metadata_key(seq_id, item.start, item.end)
                if fragment_key:
                    metadata[fragment_key] = item
            if md5:
                metadata[md5] = item
    return metadata


def te_fragment_metadata_key(seq_id: str, start: str, end: str) -> str:
    """Return the metadata key shared by extracted TE fragment BLAST subjects."""

    if not seq_id or not start or not end:
        return ""
    left = _parse_int(start)
    right = _parse_int(end)
    if left is None or right is None:
        return ""
    return "%s::te:%d-%d" % (seq_id, min(left, right), max(left, right))


def read_blast_hits(path: str | Path, te_metadata: dict[str, TEMetadata] | None = None) -> dict[str, list[TEHit]]:
    """Read BLASTN tabular output and group hits by PanPop allele MD5."""

    metadata = te_metadata or {}
    hits_by_md5: dict[str, list[TEHit]] = {}
    with Path(path).open("rt", encoding="utf-8") as handle:
        first_data_line: str | None = None
        for raw_line in handle:
            if raw_line.strip() and not raw_line.startswith("#"):
                first_data_line = raw_line.rstrip("\n")
                break
        if first_data_line is None:
            return hits_by_md5
        first_fields = first_data_line.split("\t")
        if _looks_like_blast_header(first_fields):
            rows = csv.DictReader(handle, delimiter="\t", fieldnames=[field.lower() for field in first_fields])
        else:
            rows = _iter_dict_rows([first_data_line], handle)
        for row in rows:
            qseqid = row.get("qseqid") or row.get("query") or row.get("query_id") or ""
            sseqid = row.get("sseqid") or row.get("subject") or row.get("subject_id") or ""
            allele_md5 = md5_from_identifier(qseqid)
            if not allele_md5:
                continue
            pident = _parse_float(row.get("pident") or row.get("identity"))
            qstart = _parse_int(row.get("qstart") or row.get("start"))
            qend = _parse_int(row.get("qend") or row.get("end"))
            if pident is None or qstart is None or qend is None:
                continue
            te_metadata_item = te_metadata_for_subject(sseqid, metadata)
            hits_by_md5.setdefault(allele_md5, []).append(
                TEHit(
                    start=qstart,
                    end=qend,
                    identity=pident,
                    subject_id=sseqid,
                    family="" if te_metadata_item is None else te_metadata_item.family,
                    te_class="" if te_metadata_item is None else te_metadata_item.te_class,
                    metadata=te_metadata_item,
                )
            )
    return hits_by_md5


def te_metadata_for_subject(subject_id: str, metadata: dict[str, TEMetadata]) -> TEMetadata | None:
    """Return TE metadata for a BLAST subject ID or extracted TE fragment ID."""

    if subject_id in metadata:
        return metadata[subject_id]
    if "::te:" in subject_id:
        fragment_key = subject_id.rsplit(":", 1)[0]
        if fragment_key in metadata:
            return metadata[fragment_key]
        source_id = subject_id.split("::te:", 1)[0]
        if source_id in metadata:
            return metadata[source_id]
        source_md5 = md5_from_identifier(source_id)
        if source_md5 in metadata:
            return metadata[source_md5]
    subject_md5 = md5_from_identifier(subject_id)
    if subject_md5 in metadata:
        return metadata[subject_md5]
    return None


def md5_from_identifier(identifier: str) -> str:
    fields = identifier.split("|")
    if len(fields) >= 8 and fields[7]:
        return fields[7]
    if re.fullmatch(r"[0-9a-fA-F]{32}", identifier):
        return identifier.lower()
    return ""


def classify_allele(
    sequence: str,
    hits: Sequence[TEHit],
    *,
    md5: str | None = None,
    allele_length: int | None = None,
    min_te_coverage: float = 0.60,
    min_identity: float = 80.0,
    min_te_covered_bp: int = 40,
) -> AlleleEvidence:
    allele_length = len(sequence) if allele_length is None else allele_length
    md5 = sequence_md5(sequence) if md5 is None else md5
    if allele_length == 0 or not hits:
        return AlleleEvidence(md5, allele_length, 0, 0.0, None, False)

    clipped_intervals: list[tuple[int, int]] = []
    weighted_identity_sum = 0.0
    aligned_bp_sum = 0
    family_bp: dict[tuple[str, str], int] = {}
    supporting_te_annotations: list[str] = []
    seen_supporting_annotations: set[str] = set()
    for hit in hits:
        start = max(min(hit.start, hit.end), 1)
        end = min(max(hit.start, hit.end), allele_length)
        if end < start:
            continue
        length = end - start + 1
        clipped_intervals.append((start, end))
        weighted_identity_sum += length * hit.identity
        aligned_bp_sum += length
        key = (hit.family, hit.te_class)
        family_bp[key] = family_bp.get(key, 0) + length
        supporting_annotation = format_supporting_te_annotation(hit, start, end)
        if supporting_annotation and supporting_annotation not in seen_supporting_annotations:
            supporting_te_annotations.append(supporting_annotation)
            seen_supporting_annotations.add(supporting_annotation)

    te_covered_bp = union_interval_length(clipped_intervals)
    coverage = te_covered_bp / allele_length if allele_length else 0.0
    weighted_identity = weighted_identity_sum / aligned_bp_sum if aligned_bp_sum else None
    is_te = (
        coverage >= min_te_coverage
        and weighted_identity is not None
        and weighted_identity >= min_identity
        and te_covered_bp >= min_te_covered_bp
    )
    family, te_class = ("", "")
    if family_bp:
        family, te_class = max(family_bp.items(), key=lambda item: item[1])[0]
    return AlleleEvidence(
        md5=md5,
        allele_length=allele_length,
        te_covered_bp=te_covered_bp,
        coverage=coverage,
        weighted_identity=weighted_identity,
        is_te=is_te,
        family=family,
        te_class=te_class,
        supporting_te_annotations=";".join(supporting_te_annotations),
    )


def format_supporting_te_annotation(hit: TEHit, query_start: int, query_end: int) -> str:
    """Format one supporting TE annotation as a compact report field value."""

    fields = [
        "subject=%s" % sanitize_report_value(hit.subject_id),
        "q=%d-%d" % (query_start, query_end),
        "identity=%.6g" % hit.identity,
    ]
    metadata = hit.metadata
    if metadata is not None:
        if metadata.seq_id:
            fields.append("seq_id=%s" % sanitize_report_value(metadata.seq_id))
        if metadata.start or metadata.end:
            fields.append("te=%s-%s" % (sanitize_report_value(metadata.start), sanitize_report_value(metadata.end)))
        if metadata.strand:
            fields.append("strand=%s" % sanitize_report_value(metadata.strand))
        if metadata.feature_type:
            fields.append("type=%s" % sanitize_report_value(metadata.feature_type))
        if metadata.family:
            fields.append("family=%s" % sanitize_report_value(metadata.family))
        if metadata.te_class:
            fields.append("class=%s" % sanitize_report_value(metadata.te_class))
        if metadata.attributes:
            fields.append("attributes=%s" % sanitize_report_value(metadata.attributes))
    else:
        if hit.family:
            fields.append("family=%s" % sanitize_report_value(hit.family))
        if hit.te_class:
            fields.append("class=%s" % sanitize_report_value(hit.te_class))
    return ",".join(fields)


def sanitize_report_value(value: str) -> str:
    return value.replace("\t", " ").replace("\n", " ").replace("\r", " ").replace(";", "%3B")


def union_interval_length(intervals: Iterable[tuple[int, int]]) -> int:
    sorted_intervals = sorted(intervals)
    if not sorted_intervals:
        return 0
    total = 0
    current_start, current_end = sorted_intervals[0]
    for start, end in sorted_intervals[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
            continue
        total += current_end - current_start + 1
        current_start, current_end = start, end
    total += current_end - current_start + 1
    return total


def iter_tip_vcf_lines(
    panpop_vcf: str | Path,
    hits_by_md5: dict[str, list[TEHit]],
    *,
    allele_index: dict[tuple[str, int, int, str, str, int], PanpopAlleleIndexEntry] | None = None,
    min_te_coverage: float = 0.60,
    min_identity: float = 80.0,
    min_te_covered_bp: int = 40,
) -> Iterator[ProcessedVcfRecord]:
    with Path(panpop_vcf).open("rt", encoding="utf-8") as handle:
        inserted_info = False
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("##"):
                yield ProcessedVcfRecord(line=line)
                continue
            if line.startswith("#CHROM"):
                for info_line in TIP_INFO_LINES:
                    yield ProcessedVcfRecord(line=info_line)
                inserted_info = True
                yield ProcessedVcfRecord(line=line)
                continue
            if line.startswith("#"):
                yield ProcessedVcfRecord(line=line)
                continue
            if not inserted_info:
                for info_line in TIP_INFO_LINES:
                    yield ProcessedVcfRecord(line=info_line)
                inserted_info = True
            yield process_vcf_record_line(
                line,
                hits_by_md5,
                allele_index=allele_index,
                min_te_coverage=min_te_coverage,
                min_identity=min_identity,
                min_te_covered_bp=min_te_covered_bp,
            )


def process_vcf_record_line(
    line: str,
    hits_by_md5: dict[str, list[TEHit]],
    *,
    allele_index: dict[tuple[str, int, int, str, str, int], PanpopAlleleIndexEntry] | None = None,
    min_te_coverage: float = 0.60,
    min_identity: float = 80.0,
    min_te_covered_bp: int = 40,
) -> ProcessedVcfRecord:
    fields = line.split("\t")
    if len(fields) < 8:
        return ProcessedVcfRecord(line=None)
    chrom, pos, record_id, ref, alt_raw = fields[:5]
    alts = alt_raw.split(",") if alt_raw else []
    records = parse_vcf_record_line(line, genome="panpop")

    alt_evidence: dict[int, AlleleEvidence] = {}
    ref_evidence: AlleleEvidence | None = None
    ref_is_te = False
    te_alt_numbers: list[int] = []
    retained_alt_numbers: list[int] = []
    reports: list[AlleleReportRow] = []

    for record in records:
        old_alt_number = record.allele_index + 1
        if record.svtype == "INS":
            index_entry = lookup_panpop_allele_index_entry(allele_index, record.chrom, record.pos, record.end, record.svtype, "alt", old_alt_number)
            md5 = sequence_md5(record.alt) if index_entry is None else index_entry.md5
            evidence = classify_allele(
                record.alt,
                hits_by_md5.get(md5, []),
                md5=md5,
                allele_length=None if index_entry is None else index_entry.length,
                min_te_coverage=min_te_coverage,
                min_identity=min_identity,
                min_te_covered_bp=min_te_covered_bp,
            )
            alt_evidence[old_alt_number] = evidence
            if evidence.is_te:
                te_alt_numbers.append(old_alt_number)
                retained_alt_numbers.append(old_alt_number)
        elif record.svtype == "DEL":
            if ref_evidence is None:
                index_entry = lookup_panpop_allele_index_entry(allele_index, record.chrom, record.pos, record.end, record.svtype, "ref", 0)
                md5 = sequence_md5(ref) if index_entry is None else index_entry.md5
                ref_evidence = classify_allele(
                    ref,
                    hits_by_md5.get(md5, []),
                    md5=md5,
                    allele_length=None if index_entry is None else index_entry.length,
                    min_te_coverage=min_te_coverage,
                    min_identity=min_identity,
                    min_te_covered_bp=min_te_covered_bp,
                )
                ref_is_te = ref_evidence.is_te
            if ref_is_te:
                retained_alt_numbers.append(old_alt_number)

    retained_alt_numbers = sorted(set(retained_alt_numbers))
    te_alt_numbers = sorted(set(te_alt_numbers))
    if not retained_alt_numbers:
        return ProcessedVcfRecord(line=None, reports=reports)

    old_to_new = {old: new for new, old in enumerate(retained_alt_numbers, start=1)}
    new_alts = [alts[old - 1] for old in retained_alt_numbers if 0 <= old - 1 < len(alts)]
    if not new_alts:
        return ProcessedVcfRecord(line=None, reports=reports)

    if ref_evidence is not None:
        reports.append(
            build_report_row(
                record_id=record_id,
                chrom=chrom,
                pos=pos,
                svtype="DEL",
                allele_role="ref",
                original_allele_number=0,
                new_allele_number="0",
                evidence=ref_evidence,
                retained=True,
            )
        )
    for old_alt_number, evidence in sorted(alt_evidence.items()):
        reports.append(
            build_report_row(
                record_id=record_id,
                chrom=chrom,
                pos=pos,
                svtype="INS",
                allele_role="alt",
                original_allele_number=old_alt_number,
                new_allele_number=str(old_to_new[old_alt_number]) if old_alt_number in old_to_new else ".",
                evidence=evidence,
                retained=old_alt_number in old_to_new,
            )
        )

    fields[4] = ",".join(new_alts)
    fields[7] = append_tip_info(fields[7], ref_is_te=ref_is_te, te_alt_numbers=te_alt_numbers, retained_alt_numbers=retained_alt_numbers)
    if len(fields) > 9:
        fields[9:] = [rewrite_sample_field(fields[8], sample_field, old_to_new) for sample_field in fields[9:]]
    return ProcessedVcfRecord(line="\t".join(fields), reports=reports)


def lookup_panpop_allele_index_entry(
    allele_index: dict[tuple[str, int, int, str, str, int], PanpopAlleleIndexEntry] | None,
    chrom: str,
    pos: int,
    end: int,
    svtype: str,
    allele_role: str,
    allele_number: int,
) -> PanpopAlleleIndexEntry | None:
    if allele_index is None:
        return None
    return allele_index.get(panpop_allele_index_key(chrom, pos, end, svtype, allele_role, allele_number))


def build_report_row(
    *,
    record_id: str,
    chrom: str,
    pos: str,
    svtype: str,
    allele_role: str,
    original_allele_number: int,
    new_allele_number: str,
    evidence: AlleleEvidence,
    retained: bool,
) -> AlleleReportRow:
    return AlleleReportRow(
        record_id=record_id,
        chrom=chrom,
        pos=pos,
        svtype=svtype,
        allele_role=allele_role,
        original_allele_number=original_allele_number,
        new_allele_number=new_allele_number,
        allele_length=evidence.allele_length,
        md5=evidence.md5,
        te_covered_bp=evidence.te_covered_bp,
        te_coverage=evidence.coverage,
        weighted_identity="" if evidence.weighted_identity is None else "%.6g" % evidence.weighted_identity,
        is_te=evidence.is_te,
        retained=retained,
        family=evidence.family,
        te_class=evidence.te_class,
        supporting_te_annotations=evidence.supporting_te_annotations,
    )


def append_tip_info(info: str, *, ref_is_te: bool, te_alt_numbers: Sequence[int], retained_alt_numbers: Sequence[int]) -> str:
    items = [] if info in {"", "."} else [info]
    items.append("TIPMAP")
    items.append("TIP_TE_REF=%d" % (1 if ref_is_te else 0))
    items.append("TIP_TE_ALTS=%s" % (",".join(str(value) for value in te_alt_numbers) if te_alt_numbers else "."))
    items.append("TIP_RETAINED_ALTS=%s" % ",".join(str(value) for value in retained_alt_numbers))
    return ";".join(items)


def rewrite_sample_field(format_field: str, sample_field: str, old_to_new: dict[int, int]) -> str:
    format_keys = format_field.split(":")
    sample_values = sample_field.split(":")
    try:
        gt_index = format_keys.index("GT")
    except ValueError:
        return sample_field
    if gt_index >= len(sample_values):
        return sample_field
    sample_values[gt_index] = rewrite_genotype(sample_values[gt_index], old_to_new)
    return ":".join(sample_values)


def rewrite_genotype(genotype: str, old_to_new: dict[int, int]) -> str:
    delimiter = "|" if "|" in genotype else "/"
    alleles = re.split(r"[/|]", genotype)
    rewritten = []
    for allele in alleles:
        if allele in {".", ""}:
            rewritten.append(allele or ".")
            continue
        try:
            allele_number = int(allele)
        except ValueError:
            rewritten.append(".")
            continue
        if allele_number == 0:
            rewritten.append("0")
        else:
            rewritten.append(str(old_to_new.get(allele_number, 0)))
    return delimiter.join(rewritten)


def write_tip_outputs(
    *,
    panpop_vcf: str | Path,
    blast_tsv: str | Path,
    output_vcf: str | Path,
    te_annotations: str | Path | None = None,
    allele_report: str | Path | None = None,
    allele_index: dict[tuple[str, int, int, str, str, int], PanpopAlleleIndexEntry] | None = None,
    min_te_coverage: float = 0.60,
    min_identity: float = 80.0,
    min_te_covered_bp: int = 40,
) -> tuple[int, int]:
    metadata = read_te_metadata(te_annotations)
    hits_by_md5 = read_blast_hits(blast_tsv, metadata)
    output_path = Path(output_vcf)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_rows: list[AlleleReportRow] = []
    record_count = 0
    with output_path.open("wt", encoding="utf-8", newline="\n") as handle:
        for processed in iter_tip_vcf_lines(
            panpop_vcf,
            hits_by_md5,
            allele_index=allele_index,
            min_te_coverage=min_te_coverage,
            min_identity=min_identity,
            min_te_covered_bp=min_te_covered_bp,
        ):
            report_rows.extend(processed.reports)
            if processed.line is None:
                continue
            handle.write(processed.line + "\n")
            if not processed.line.startswith("#"):
                record_count += 1
    if allele_report is not None:
        write_allele_report(report_rows, allele_report)
    return record_count, len(report_rows)


def run_classification_workflow(
    *,
    panpop_vcf: str | Path,
    output_vcf: str | Path,
    sv_fasta: str | Path | None = None,
    te_annotations: str | Path | None = None,
    blast_tsv: str | Path | None = None,
    allele_report: str | Path | None = None,
    workdir: str | Path | None = None,
    panpop_alleles_fasta: str | Path | None = None,
    te_fragments_fasta: str | Path | None = None,
    makeblastdb: str = "makeblastdb",
    blast_db_prefix: str | Path | None = None,
    blastn: str = "blastn",
    blast_threads: int = 1,
    blast_evalue: float = 1e-5,
    blast_args: Sequence[str] = (),
    min_panpop_allele_length: int = 1,
    min_te_fragment_length: int = 1,
    fasta_line_width: int = 80,
    min_te_coverage: float = 0.60,
    min_identity: float = 80.0,
    min_te_covered_bp: int = 40,
    runner: CommandRunner | None = None,
) -> tuple[int, int]:
    output_path = Path(output_vcf)
    work_path = Path(workdir) if workdir is not None else output_path.with_suffix(".tipmap_work")
    work_path.mkdir(parents=True, exist_ok=True)
    blast_path = Path(blast_tsv) if blast_tsv is not None else work_path / "panpop_allele_vs_te.tsv"
    allele_index: dict[tuple[str, int, int, str, str, int], PanpopAlleleIndexEntry] | None = None

    if blast_tsv is None:
        if sv_fasta is None or te_annotations is None:
            msg = "--sv-fasta and --te-annotations are required when --blast-tsv is not provided"
            raise ValueError(msg)
        allele_path = Path(panpop_alleles_fasta) if panpop_alleles_fasta is not None else work_path / "panpop_alleles.fa"
        allele_index_path = work_path / "panpop_alleles.index.tsv"
        fragment_path = Path(te_fragments_fasta) if te_fragments_fasta is not None else work_path / "te_fragments.fa"
        dedup_fragment_path = work_path / "te_fragments.dedup.fa"
        dedup_metadata_path = work_path / "te_fragments.dedup.metadata.tsv"
        db_prefix = Path(blast_db_prefix) if blast_db_prefix is not None else work_path / "te_fragments_db"
        allele_index = write_panpop_te_alleles_and_index(
            panpop_vcf,
            allele_path,
            allele_index_path,
            min_length=min_panpop_allele_length,
            line_width=fasta_line_width,
        )
        write_te_fragments(
            sv_fasta,
            te_annotations,
            fragment_path,
            min_length=min_te_fragment_length,
            line_width=fasta_line_width,
        )
        write_deduplicated_te_fragments(
            fragment_path,
            dedup_fragment_path,
            dedup_metadata_path,
            line_width=fasta_line_width,
        )
        run_makeblastdb(
            fasta=dedup_fragment_path,
            db_prefix=db_prefix,
            makeblastdb=makeblastdb,
            runner=runner,
        )
        run_blastn(
            query_fasta=allele_path,
            blast_db=db_prefix,
            output=blast_path,
            blastn=blastn,
            blast_threads=blast_threads,
            blast_evalue=blast_evalue,
            extra_args=blast_args,
            runner=runner,
        )

    return write_tip_outputs(
        panpop_vcf=panpop_vcf,
        blast_tsv=blast_path,
        te_annotations=te_annotations,
        output_vcf=output_vcf,
        allele_report=allele_report,
        allele_index=allele_index,
        min_te_coverage=min_te_coverage,
        min_identity=min_identity,
        min_te_covered_bp=min_te_covered_bp,
    )


def write_allele_report(rows: Sequence[AlleleReportRow], output: str | Path) -> None:
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["class" if field == "te_class" else field for field in AlleleReportRow.__dataclass_fields__]
    with output_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            raw_row = row.__dict__.copy()
            raw_row["class"] = raw_row.pop("te_class")
            writer.writerow(raw_row)


def _looks_like_blast_header(fields: Sequence[str]) -> bool:
    lower = {field.lower() for field in fields}
    return "qseqid" in lower and ("pident" in lower or "identity" in lower)


def _iter_dict_rows(first_lines: Sequence[str], handle: Iterable[str]) -> Iterator[dict[str, str]]:
    for line in list(first_lines) + [raw.rstrip("\n") for raw in handle if raw.strip() and not raw.startswith("#")]:
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        yield {column: fields[index] for index, column in enumerate(BLAST_OUTFMT6_COLUMNS) if index < len(fields)}


def _parse_int(raw_value: str | None) -> int | None:
    if raw_value in {None, "", "."}:
        return None
    return int(float(raw_value))


def _parse_float(raw_value: str | None) -> float | None:
    if raw_value in {None, "", "."}:
        return None
    value = float(raw_value)
    return value * 100.0 if value <= 1.0 else value


def _default_runner(command: Sequence[str], cwd: Path | None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(command),
        cwd=None if cwd is None else str(cwd),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()
    try:
        run_classification_workflow(
            panpop_vcf=args.panpop_vcf,
            sv_fasta=args.sv_fasta,
            te_annotations=args.te_annotations,
            blast_tsv=args.blast_tsv,
            output_vcf=args.output_vcf,
            allele_report=args.allele_report,
            workdir=args.workdir,
            panpop_alleles_fasta=args.panpop_alleles_fasta,
            te_fragments_fasta=args.te_fragments_fasta,
            makeblastdb=args.makeblastdb,
            blast_db_prefix=args.blast_db_prefix,
            blastn=args.blastn,
            blast_threads=args.blast_threads,
            blast_evalue=args.blast_evalue,
            blast_args=args.blast_arg,
            min_panpop_allele_length=args.min_panpop_allele_length,
            min_te_fragment_length=args.min_te_fragment_length,
            fasta_line_width=args.fasta_line_width,
            min_te_coverage=args.min_te_coverage,
            min_identity=args.min_identity,
            min_te_covered_bp=args.min_te_covered_bp,
        )
    except ValueError as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())








