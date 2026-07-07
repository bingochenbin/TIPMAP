"""VCF parsing utilities for PanPop SV-VCF and pairwise genome-to-reference VCF records."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
import gzip
from pathlib import Path
import re
from typing import Dict, List, Optional, TextIO

from tipmap.lib.models import AltKind, SVRecord, VCFSource

SUPPORTED_SVTYPES = {"INS", "DEL"}
BREAKEND_MARKERS = {"[", "]"}
TRUE_FLAG = "True"


class VCFParseError(ValueError):
    """Raised when a VCF line cannot be parsed into a TIPMap SV record."""


def parse_vcf(
    path: str | Path,
    *,
    genome: Optional[str] = None,
    vcf_source: VCFSource = "unknown",
    use_header_sample: bool = True,
    parse_sample_genotypes: bool = False,
) -> Iterator[SVRecord]:
    """Yield TIPMap SV records from a PanPop or pairwise VCF file.

    The parser reads plain text ``.vcf`` files and gzip-compressed ``.vcf.gz``
    files. If ``genome`` is not provided, pairwise VCF parsing can use the first
    sample name from the VCF header; PanPop parsing uses the file stem by default.
    """

    vcf_path = Path(path)
    with _open_text(vcf_path) as handle:
        yield from parse_vcf_lines(
            handle,
            genome=genome,
            fallback_genome=_vcf_stem(vcf_path),
            vcf_source=vcf_source,
            use_header_sample=use_header_sample,
            parse_sample_genotypes=parse_sample_genotypes,
        )


def parse_panpop_vcf(path: str | Path, *, genome: Optional[str] = None) -> Iterator[SVRecord]:
    """Yield records from a PanPop-derived merged SV-VCF.

    PanPop records define merged SV coordinates, IDs, and sample genotypes after
    left alignment and normalization. Their REF/ALT columns can be used to infer
    INS/DEL type and length, but pairwise VCFs remain the preferred source for
    true ALT sequences.
    """

    yield from parse_vcf(
        path,
        genome=genome,
        vcf_source="panpop",
        use_header_sample=False,
        parse_sample_genotypes=True,
    )


def parse_pairwise_vcf(path: str | Path, *, genome: Optional[str] = None) -> Iterator[SVRecord]:
    """Yield records from a pairwise genome-to-reference VCF.

    Pairwise alignment VCFs are treated as the preferred source for true ALT
    sequence whenever the ALT column contains a nucleotide sequence. ``ASM_*``
    INFO fields, when present, describe ALT coordinates in the sample assembly
    and are preserved in ``SVRecord.info``.
    """

    yield from parse_vcf(path, genome=genome, vcf_source="pairwise")


def parse_vcf_lines(
    lines: Iterable[str],
    *,
    genome: Optional[str] = None,
    fallback_genome: str = "unknown",
    vcf_source: VCFSource = "unknown",
    use_header_sample: bool = True,
    parse_sample_genotypes: bool = False,
) -> Iterator[SVRecord]:
    """Yield parsed VCF records from PanPop or pairwise VCF text lines."""

    active_genome = genome
    sample_names: List[str] = []
    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            sample_names = _sample_names_from_header(line)
            if active_genome is None and use_header_sample:
                active_genome = _sample_name_from_header(line)
            continue
        if line.startswith("#"):
            continue

        record_genome = active_genome or fallback_genome
        try:
            yield from parse_vcf_record_line(
                line,
                genome=record_genome,
                vcf_source=vcf_source,
                sample_names=sample_names,
                parse_sample_genotypes=parse_sample_genotypes,
            )
        except Exception as exc:
            msg = "Failed to parse VCF line %d: %s" % (line_number, exc)
            raise VCFParseError(msg) from exc


def parse_vcf_record_line(
    line: str,
    *,
    genome: str,
    vcf_source: VCFSource = "unknown",
    sample_names: Optional[List[str]] = None,
    parse_sample_genotypes: bool = False,
) -> List[SVRecord]:
    """Parse one non-header VCF record line.

    Multi-allelic records are expanded into one :class:`SVRecord` per ALT
    allele. Per-allele ``SVTYPE`` and ``SVLEN`` values are selected by ALT index
    when those INFO fields contain comma-separated values.
    """

    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8:
        msg = "VCF records must contain at least 8 tab-delimited columns"
        raise VCFParseError(msg)

    chrom, pos_raw, source_id_raw, ref, alt_raw, qual_raw, filter_raw, info_raw = fields[:8]
    pos = _parse_required_int(pos_raw, "POS")
    info = parse_info_field(info_raw)
    alts = alt_raw.split(",")
    sample_genotypes = (
        _parse_sample_genotypes(fields, sample_names or []) if parse_sample_genotypes else {}
    )

    records = []
    for allele_index, alt in enumerate(alts):
        alt_kind = detect_alt_kind(alt)
        svtype = _detect_svtype(info, ref, alt, allele_index)
        end = _detect_end(info, pos, ref, alt, svtype, allele_index)
        svlen = _detect_svlen(info, pos, end, ref, alt, svtype, allele_index)
        records.append(
            SVRecord(
                genome=genome,
                chrom=chrom,
                pos=pos,
                end=end,
                svtype=svtype,
                svlen=svlen,
                ref=ref,
                alt=alt,
                source_id=None if source_id_raw == "." else source_id_raw,
                qual=_parse_optional_float(qual_raw),
                filter=None if filter_raw in {".", ""} else filter_raw,
                info=info,
                sample_genotypes=sample_genotypes,
                allele_index=allele_index,
                vcf_source=vcf_source,
                alt_kind=alt_kind,
            )
        )
    return records


def parse_info_field(raw_info: str) -> Dict[str, str]:
    """Parse a VCF INFO field into a string dictionary."""

    if raw_info in {"", "."}:
        return {}

    info: Dict[str, str] = {}
    for item in raw_info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = TRUE_FLAG
    return info


def detect_alt_kind(alt: str) -> AltKind:
    """Classify a VCF ALT allele by whether it contains true sequence."""

    if alt in {"", "."}:
        return "missing"
    if alt.startswith("<") and alt.endswith(">"):
        return "symbolic"
    if any(marker in alt for marker in BREAKEND_MARKERS):
        return "breakend"
    return "sequence"


def _detect_svtype(
    info: Dict[str, str],
    ref: str,
    alt: str,
    allele_index: int,
) -> str:
    info_svtype = _select_info_value(info, "SVTYPE", allele_index)
    if info_svtype:
        return _normalize_svtype(info_svtype)

    symbolic_svtype = _svtype_from_symbolic_alt(alt)
    if symbolic_svtype:
        return symbolic_svtype

    if any(marker in alt for marker in BREAKEND_MARKERS):
        return "BND"

    length_delta = len(alt) - len(ref)
    if length_delta > 0:
        return "INS"
    if length_delta < 0:
        return "DEL"
    return "UNK"


def _detect_end(
    info: Dict[str, str],
    pos: int,
    ref: str,
    alt: str,
    svtype: str,
    allele_index: int,
) -> int:
    info_end = _parse_optional_int(_select_info_value(info, "END", allele_index))
    if info_end is not None:
        return info_end

    info_svlen = _parse_optional_int(_select_info_value(info, "SVLEN", allele_index))
    if info_svlen is not None and svtype == "DEL":
        return pos + abs(info_svlen) - 1

    if _is_sequence_alt(alt) and len(ref) > 1:
        return pos + len(ref) - 1

    return pos


def _detect_svlen(
    info: Dict[str, str],
    pos: int,
    end: int,
    ref: str,
    alt: str,
    svtype: str,
    allele_index: int,
) -> int:
    info_svlen = _parse_optional_int(_select_info_value(info, "SVLEN", allele_index))
    if info_svlen is not None:
        return info_svlen

    if _is_sequence_alt(alt):
        delta = len(alt) - len(ref)
        if delta != 0:
            return delta

    span = max(end - pos + 1, 0)
    if svtype == "DEL":
        return -span
    if svtype == "INS":
        return span
    return 0


def _normalize_svtype(raw_svtype: str) -> str:
    cleaned = raw_svtype.strip().upper()
    cleaned = cleaned.strip("<>")
    cleaned = cleaned.split(":", 1)[0]
    cleaned = re.sub(r"[^A-Z0-9_]+", "", cleaned)
    if cleaned in SUPPORTED_SVTYPES:
        return cleaned
    return "UNK"


def _svtype_from_symbolic_alt(alt: str) -> Optional[str]:
    if not (alt.startswith("<") and alt.endswith(">")):
        return None
    return _normalize_svtype(alt)


def _is_sequence_alt(alt: str) -> bool:
    return detect_alt_kind(alt) == "sequence"


def _select_info_value(info: Dict[str, str], key: str, allele_index: int) -> Optional[str]:
    value = info.get(key)
    if value is None or value == ".":
        return None
    values = value.split(",")
    if allele_index < len(values):
        selected = values[allele_index]
    else:
        selected = values[0]
    return None if selected in {"", "."} else selected


def _parse_required_int(raw_value: str, field_name: str) -> int:
    value = _parse_optional_int(raw_value)
    if value is None:
        msg = "%s must be an integer" % field_name
        raise VCFParseError(msg)
    return value


def _parse_optional_int(raw_value: Optional[str]) -> Optional[int]:
    if raw_value is None or raw_value in {"", "."}:
        return None
    return int(float(raw_value))


def _parse_optional_float(raw_value: str) -> Optional[float]:
    if raw_value in {"", "."}:
        return None
    return float(raw_value)


def _sample_name_from_header(header_line: str) -> Optional[str]:
    sample_names = _sample_names_from_header(header_line)
    if sample_names:
        return sample_names[0]
    return None


def _sample_names_from_header(header_line: str) -> List[str]:
    fields = header_line.lstrip("#").split("\t")
    if len(fields) > 9:
        return fields[9:]
    return []


def _parse_sample_genotypes(fields: List[str], sample_names: List[str]) -> Dict[str, str]:
    if len(fields) < 10 or not sample_names:
        return {}

    format_keys = fields[8].split(":")
    try:
        gt_index = format_keys.index("GT")
    except ValueError:
        return {}

    genotypes: Dict[str, str] = {}
    for sample_name, sample_value in zip(sample_names, fields[9:]):
        sample_fields = sample_value.split(":")
        if gt_index < len(sample_fields):
            genotypes[sample_name] = sample_fields[gt_index]
    return genotypes


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8")
    return path.open(mode="rt", encoding="utf-8")


def _vcf_stem(path: Path) -> str:
    name = path.name
    if name.endswith(".vcf.gz"):
        return name[: -len(".vcf.gz")]
    if name.endswith(".vcf"):
        return name[: -len(".vcf")]
    return path.stem
