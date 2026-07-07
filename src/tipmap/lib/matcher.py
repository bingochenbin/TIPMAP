"""Match PanPop SV records to true ALT alleles from pairwise VCFs."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
import statistics
from typing import Dict, List, Optional, Sequence, Tuple

from tipmap.lib.fasta import IndexedFasta, extract_interval_from_source, reverse_complement
from tipmap.lib.models import ExtractedSequence, SVRecord
from tipmap.lib.utils import sequence_md5

SUPPORTED_TIP_SVTYPES = {"INS", "DEL"}


@dataclass(frozen=True)
class MatchCriteria:
    """Coordinate and length tolerances for matching pairwise VCF ALT records."""

    max_distance: int = 300
    max_length_difference: Optional[int] = None
    max_length_ratio_difference: float = 0.25


@dataclass(frozen=True)
class _IndexedPairwiseRecord:
    """A pairwise record with its original input order."""

    order: int
    record: SVRecord


@dataclass(frozen=True)
class _PairwiseRecordBucket:
    """Pairwise records for one genome/chromosome/SV type, sorted by POS."""

    positions: Tuple[int, ...]
    records: Tuple[_IndexedPairwiseRecord, ...]


class PairwiseRecordIndex:
    """Position index for fast PanPop-to-pairwise ALT lookup.

    The index stores pairwise ALT records by genome, reference chromosome, and
    inferred SV type. Within each bucket records are sorted by reference POS, so
    lookup only inspects records in the PanPop POS tolerance window before the
    full coordinate/length/genotype matching rule is applied.
    """

    def __init__(self, records: Iterable[SVRecord]) -> None:
        grouped: Dict[Tuple[str, str, str], List[_IndexedPairwiseRecord]] = {}
        locus_keys: Dict[Tuple[str, str], set[Tuple[str, str, str]]] = {}
        for order, record in enumerate(records):
            if not is_supported_tip_sv(record):
                continue
            if not record.has_real_alt_sequence:
                continue
            key = (record.genome, record.chrom, record.svtype)
            grouped.setdefault(key, []).append(_IndexedPairwiseRecord(order, record))
            locus_keys.setdefault((record.chrom, record.svtype), set()).add(key)

        self._buckets: Dict[Tuple[str, str, str], _PairwiseRecordBucket] = {}
        for key, entries in grouped.items():
            sorted_entries = tuple(sorted(entries, key=lambda entry: (entry.record.pos, entry.order)))
            self._buckets[key] = _PairwiseRecordBucket(
                positions=tuple(entry.record.pos for entry in sorted_entries),
                records=sorted_entries,
            )
        self._locus_keys: Dict[Tuple[str, str], Tuple[Tuple[str, str, str], ...]] = {
            key: tuple(sorted(keys)) for key, keys in locus_keys.items()
        }

    def find_matching_alt_records(
        self,
        panpop_record: SVRecord,
        *,
        criteria: MatchCriteria = MatchCriteria(),
    ) -> List[SVRecord]:
        """Return indexed pairwise records matching one PanPop SV."""

        if not is_supported_tip_sv(panpop_record):
            return []

        matches: List[_IndexedPairwiseRecord] = []
        for key in self._candidate_bucket_keys(panpop_record):
            bucket = self._buckets.get(key)
            if bucket is None:
                continue
            left = bisect_left(bucket.positions, panpop_record.pos - criteria.max_distance)
            right = bisect_right(bucket.positions, panpop_record.pos + criteria.max_distance)
            for entry in bucket.records[left:right]:
                if sv_records_match(panpop_record, entry.record, criteria=criteria):
                    matches.append(entry)

        matches.sort(key=lambda entry: entry.order)
        return [entry.record for entry in matches]

    def _candidate_bucket_keys(self, panpop_record: SVRecord) -> Sequence[Tuple[str, str, str]]:
        if panpop_record.sample_genotypes:
            return tuple(
                (genome, panpop_record.chrom, panpop_record.svtype)
                for genome in panpop_record.sample_genotypes
                if panpop_record.carries_alt(genome)
            )
        return self._locus_keys.get((panpop_record.chrom, panpop_record.svtype), ())


def panpop_sv_key(record: SVRecord) -> str:
    """Return a stable external key for a PanPop SV record."""

    if record.source_id:
        return record.source_id
    return "%s:%d-%d:%s:%d" % (
        record.chrom,
        record.pos,
        record.end,
        record.svtype,
        record.allele_index,
    )


def is_supported_tip_sv(record: SVRecord) -> bool:
    """Return whether the current TIPMap workflow supports this SV type."""

    return record.svtype in SUPPORTED_TIP_SVTYPES


def sv_records_match(
    panpop_record: SVRecord,
    pairwise_record: SVRecord,
    *,
    criteria: MatchCriteria = MatchCriteria(),
) -> bool:
    """Return whether a pairwise VCF record likely corresponds to a PanPop SV."""

    if not is_supported_tip_sv(panpop_record):
        return False
    if panpop_record.chrom != pairwise_record.chrom:
        return False
    if panpop_record.svtype != pairwise_record.svtype:
        return False
    if not pairwise_record.has_real_alt_sequence:
        return False
    if panpop_record.sample_genotypes and not panpop_record.carries_alt(pairwise_record.genome):
        return False
    if abs(panpop_record.pos - pairwise_record.pos) > criteria.max_distance:
        return False
    if abs(panpop_record.end - pairwise_record.end) > criteria.max_distance:
        return False

    length_difference = abs(abs(panpop_record.svlen) - abs(pairwise_record.svlen))
    if criteria.max_length_difference is not None:
        return length_difference <= criteria.max_length_difference

    denominator = max(abs(panpop_record.svlen), abs(pairwise_record.svlen), 1)
    return (length_difference / denominator) <= criteria.max_length_ratio_difference


def find_matching_alt_records(
    panpop_record: SVRecord,
    pairwise_records: Iterable[SVRecord] | PairwiseRecordIndex,
    *,
    criteria: MatchCriteria = MatchCriteria(),
) -> List[SVRecord]:
    """Return pairwise records with true ALT sequences matching one PanPop SV."""

    if isinstance(pairwise_records, PairwiseRecordIndex):
        return pairwise_records.find_matching_alt_records(panpop_record, criteria=criteria)

    return PairwiseRecordIndex(pairwise_records).find_matching_alt_records(
        panpop_record,
        criteria=criteria,
    )


def extract_ref_sequence(
    panpop_record: SVRecord,
    reference_sequences: Mapping[str, str] | IndexedFasta,
    *,
    flank: int = 0,
    reference_genome: str = "reference",
) -> ExtractedSequence:
    """Extract the reference-side sequence for one PanPop SV."""

    sequence = extract_interval_from_source(
        reference_sequences,
        panpop_record.chrom,
        panpop_record.pos,
        panpop_record.end,
        flank=flank,
    )
    return ExtractedSequence(
        panpop_id=panpop_sv_key(panpop_record),
        role="ref",
        genome=reference_genome,
        chrom=panpop_record.chrom,
        pos=panpop_record.pos,
        end=panpop_record.end,
        svtype=panpop_record.svtype,
        sequence=sequence,
        md5=sequence_md5(sequence),
        source_id=panpop_record.source_id,
        allele_index=panpop_record.allele_index,
    )


def should_extract_panpop_ref(record: SVRecord, min_sequence_length: int = 0) -> bool:
    """Return whether the PanPop REF sequence should be extracted."""

    if min_sequence_length < 0:
        msg = "min_sequence_length must be >= 0"
        raise ValueError(msg)
    return len(record.ref) >= min_sequence_length


def should_extract_panpop_alt(record: SVRecord, min_sequence_length: int = 0) -> bool:
    """Return whether the PanPop ALT allele should drive ALT extraction."""

    if min_sequence_length < 0:
        msg = "min_sequence_length must be >= 0"
        raise ValueError(msg)
    if record.alt_kind == "sequence":
        return len(record.alt) >= min_sequence_length
    return abs(record.svlen) >= min_sequence_length

def extract_alt_sequence(
    panpop_record: SVRecord,
    pairwise_record: SVRecord,
    *,
    query_sequences_by_genome: Mapping[str, Mapping[str, str] | IndexedFasta] | None = None,
    alt_flank: int = 0,
) -> ExtractedSequence:
    """Convert one matched pairwise VCF ALT allele to an extracted sequence."""

    sequence = pairwise_record.alt
    if alt_flank < 0:
        msg = "alt_flank must be >= 0"
        raise ValueError(msg)
    if alt_flank > 0 and query_sequences_by_genome is not None:
        query_sequence = _extract_query_flanked_alt_sequence(
            pairwise_record,
            query_sequences_by_genome,
            alt_flank=alt_flank,
        )
        if query_sequence is not None:
            sequence = query_sequence

    return ExtractedSequence(
        panpop_id=panpop_sv_key(panpop_record),
        role="alt",
        genome=pairwise_record.genome,
        chrom=pairwise_record.chrom,
        pos=pairwise_record.pos,
        end=pairwise_record.end,
        svtype=pairwise_record.svtype,
        sequence=sequence,
        md5=sequence_md5(sequence),
        source_id=pairwise_record.source_id,
        allele_index=pairwise_record.allele_index,
    )


def _extract_query_flanked_alt_sequence(
    pairwise_record: SVRecord,
    query_sequences_by_genome: Mapping[str, Mapping[str, str] | IndexedFasta],
    *,
    alt_flank: int,
) -> str | None:
    query_sequences = query_sequences_by_genome.get(pairwise_record.genome)
    if query_sequences is None:
        return None

    asm_chrom = pairwise_record.info.get("ASM_Chr")
    asm_start = _parse_positive_int(pairwise_record.info.get("ASM_Start"))
    asm_end = _parse_positive_int(pairwise_record.info.get("ASM_End"))
    asm_strand = pairwise_record.info.get("ASM_Strand", "+")
    if asm_chrom is None or asm_start is None or asm_end is None:
        return None

    query_chrom = _resolve_sequence_name(query_sequences, asm_chrom)
    if query_chrom is None:
        query_chrom = _resolve_sequence_name(query_sequences, pairwise_record.chrom)
    if query_chrom is None:
        return None

    start = min(asm_start, asm_end)
    end = max(asm_start, asm_end)
    sequence = extract_interval_from_source(query_sequences, query_chrom, start, end, flank=alt_flank)
    if asm_strand == "-":
        return reverse_complement(sequence)
    return sequence


def _parse_positive_int(value: str | None) -> int | None:
    if value is None:
        return None
    try:
        parsed = int(value)
    except ValueError:
        return None
    if parsed < 1:
        return None
    return parsed


def _resolve_sequence_name(sequences: Mapping[str, str] | IndexedFasta, chrom: str) -> str | None:
    if chrom in sequences:
        return chrom
    if chrom.startswith("Chr"):
        stripped = chrom[3:]
        if stripped in sequences:
            return stripped
    else:
        prefixed = "Chr" + chrom
        if prefixed in sequences:
            return prefixed
    return None


def iter_extracted_sequences(
    panpop_records: Iterable[SVRecord],
    pairwise_records: Iterable[SVRecord] | PairwiseRecordIndex,
    reference_sequences: Mapping[str, str] | IndexedFasta,
    *,
    criteria: MatchCriteria = MatchCriteria(),
    flank: int = 0,
    alt_flank: int = 0,
    min_panpop_sequence_length: int = 0,
    query_sequences_by_genome: Mapping[str, Mapping[str, str] | IndexedFasta] | None = None,
    reference_genome: str = "reference",
) -> Iterator[ExtractedSequence]:
    """Yield one REF plus every matched true ALT sequence for PanPop INS/DEL records."""

    if isinstance(pairwise_records, PairwiseRecordIndex):
        pairwise_index = pairwise_records
    else:
        pairwise_index = PairwiseRecordIndex(pairwise_records)
    if min_panpop_sequence_length < 0:
        msg = "min_panpop_sequence_length must be >= 0"
        raise ValueError(msg)
    for panpop_record in panpop_records:
        if not is_supported_tip_sv(panpop_record):
            continue
        if should_extract_panpop_ref(panpop_record, min_panpop_sequence_length):
            yield extract_ref_sequence(
                panpop_record,
                reference_sequences,
                flank=flank,
                reference_genome=reference_genome,
            )
        if not should_extract_panpop_alt(panpop_record, min_panpop_sequence_length):
            continue
        for pairwise_record in find_matching_alt_records(
            panpop_record,
            pairwise_index,
            criteria=criteria,
        ):
            yield extract_alt_sequence(
                panpop_record,
                pairwise_record,
                query_sequences_by_genome=query_sequences_by_genome,
                alt_flank=alt_flank,
            )


def choose_representative_alt(records: Iterable[ExtractedSequence]) -> Optional[ExtractedSequence]:
    """Choose the true ALT whose length is closest to the median ALT length."""

    alt_records = [record for record in records if record.role == "alt"]
    if not alt_records:
        return None
    median_length = statistics.median(record.length for record in alt_records)
    return min(
        alt_records,
        key=lambda record: (
            abs(record.length - median_length),
            record.length,
            record.genome,
            record.md5,
        ),
    )


__all__ = [
    "SUPPORTED_TIP_SVTYPES",
    "MatchCriteria",
    "PairwiseRecordIndex",
    "panpop_sv_key",
    "is_supported_tip_sv",
    "sv_records_match",
    "find_matching_alt_records",
    "extract_ref_sequence",
    "should_extract_panpop_ref",
    "should_extract_panpop_alt",
    "extract_alt_sequence",
    "iter_extracted_sequences",
    "choose_representative_alt",
]


