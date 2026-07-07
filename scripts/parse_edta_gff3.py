"""Parse EDTA GFF3 outputs into a TIPMap TE annotation TSV."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable, Iterator, Sequence

GFF3_PATTERNS = ("*EDTA*.gff3", "*TEanno*.gff3", "*.gff3")
TSV_HEADER = "seq_id\tmd5\tchrom\tstart\tend\tstrand\tsource\ttype\tfamily\tsuperfamily\tattributes\n"


@dataclass(frozen=True)
class TEAnnotationRow:
    """One parsed EDTA annotation row for TSV export."""

    seq_id: str
    md5: str
    chrom: str
    start: int
    end: int
    strand: str
    source: str
    feature_type: str
    family: str
    superfamily: str
    attributes: str


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gff3",
        action="append",
        default=[],
        help="EDTA GFF3 file to parse. Repeat for multiple files.",
    )
    parser.add_argument(
        "--edta-dir",
        action="append",
        default=[],
        help="EDTA output directory to scan recursively for GFF3 files. Repeat for multiple directories.",
    )
    parser.add_argument("--output", required=True, help="Output normalized TE annotation TSV path.")
    return parser


def collect_gff3_paths(gff3_files: Sequence[str | Path], edta_dirs: Sequence[str | Path]) -> list[Path]:
    """Collect explicit and directory-discovered GFF3 paths in deterministic order."""

    seen: set[Path] = set()
    paths: list[Path] = []
    for file_path in gff3_files:
        path = Path(file_path)
        if not path.is_file():
            msg = "GFF3 file not found: %s" % path
            raise FileNotFoundError(msg)
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            paths.append(path)

    for directory in edta_dirs:
        dir_path = Path(directory)
        if not dir_path.is_dir():
            msg = "EDTA directory not found: %s" % dir_path
            raise FileNotFoundError(msg)
        for pattern in GFF3_PATTERNS:
            for path in sorted(dir_path.rglob(pattern)):
                if not path.is_file():
                    continue
                resolved = path.resolve()
                if resolved in seen:
                    continue
                seen.add(resolved)
                paths.append(path)

    return sorted(paths, key=lambda path: str(path))


def parse_edta_gff3(path: str | Path) -> Iterator[TEAnnotationRow]:
    """Parse one EDTA GFF3 file into normalized annotation rows."""

    gff_path = Path(path)
    with gff_path.open("rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                continue
            seq_id, source, feature_type, start, end, _score, strand, _phase, raw_attributes = fields
            attributes = parse_gff3_attributes(raw_attributes)
            family, superfamily = classify_te_attributes(attributes, feature_type)
            yield TEAnnotationRow(
                seq_id=seq_id,
                md5=md5_from_tipmap_header(seq_id),
                chrom=chromosome_from_tipmap_header(seq_id),
                start=int(start),
                end=int(end),
                strand=strand,
                source=source,
                feature_type=feature_type,
                family=family,
                superfamily=superfamily,
                attributes=raw_attributes,
            )


def parse_gff3_attributes(raw_attributes: str) -> dict[str, str]:
    """Parse a GFF3 attribute field into a dictionary."""

    attributes: dict[str, str] = {}
    for item in raw_attributes.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes


def classify_te_attributes(attributes: dict[str, str], fallback: str) -> tuple[str, str]:
    """Extract TE family and superfamily from common EDTA attribute keys."""

    classification = (
        attributes.get("Classification")
        or attributes.get("classification")
        or attributes.get("Class")
        or attributes.get("Name")
        or attributes.get("ID")
        or fallback
    )
    parts = re.split(r"[/|:]", classification)
    cleaned = [part for part in parts if part]
    if len(cleaned) >= 2:
        return cleaned[-1], cleaned[0]
    return classification, ""


def chromosome_from_tipmap_header(header_name: str) -> str:
    """Extract chromosome from a TIPMap FASTA-derived GFF3 sequence ID."""

    fields = header_name.split("|")
    if len(fields) >= 4 and fields[3]:
        return fields[3]
    return "unknown"


def md5_from_tipmap_header(header_name: str) -> str:
    """Extract sequence MD5 from a TIPMap FASTA-derived GFF3 sequence ID."""

    fields = header_name.split("|")
    if len(fields) >= 8:
        return fields[7]
    return ""


def iter_annotation_rows(paths: Iterable[str | Path]) -> Iterator[TEAnnotationRow]:
    """Yield parsed annotation rows from all GFF3 paths."""

    for path in paths:
        yield from parse_edta_gff3(path)


def write_annotation_tsv(rows: Iterable[TEAnnotationRow], output: str | Path) -> int:
    """Write normalized TE annotation rows to TSV and return row count."""

    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with output_path.open("wt", encoding="utf-8", newline="\n") as handle:
        handle.write(TSV_HEADER)
        for row in rows:
            handle.write(
                "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n"
                % (
                    row.seq_id,
                    row.md5,
                    row.chrom,
                    row.start,
                    row.end,
                    row.strand,
                    row.source,
                    row.feature_type,
                    row.family,
                    row.superfamily,
                    row.attributes,
                )
            )
            count += 1
    return count


def run_parse_workflow(
    *,
    gff3_files: Sequence[str | Path] = (),
    edta_dirs: Sequence[str | Path] = (),
    output: str | Path,
) -> int:
    """Collect EDTA GFF3 files, parse them, and write a normalized TSV."""

    paths = collect_gff3_paths(gff3_files, edta_dirs)
    if not paths:
        msg = "No EDTA GFF3 files found. Provide --gff3 or --edta-dir with matching files."
        raise FileNotFoundError(msg)
    return write_annotation_tsv(iter_annotation_rows(paths), output)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    if not args.gff3 and not args.edta_dir:
        parser.error("at least one --gff3 or --edta-dir is required")
    try:
        run_parse_workflow(gff3_files=args.gff3, edta_dirs=args.edta_dir, output=args.output)
    except FileNotFoundError as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
