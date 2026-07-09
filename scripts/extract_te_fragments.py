"""Extract TE fragment sequences from TIPMap FASTA using TE annotation TSV coordinates."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterator

from tipmap.lib.fasta import FastaRecord, read_fasta, reverse_complement, write_fasta_records
from tipmap.lib.utils import sequence_md5


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sv-fasta", required=True, help="FASTA that was annotated by EDTA/TE tools.")
    parser.add_argument("--te-annotations", required=True, help="TE annotation TSV from annotate_te.py or parse_edta_gff3.py.")
    parser.add_argument("--output", required=True, help="Output TE fragment FASTA path.")
    parser.add_argument("--min-length", type=int, default=1, help="Minimum fragment length to emit. Default: 1.")
    parser.add_argument("--line-width", type=int, default=80, help="Output FASTA line width. Default: 80.")
    return parser


def iter_te_fragments(
    sv_fasta: str | Path,
    te_annotations: str | Path,
    *,
    min_length: int = 1,
) -> Iterator[FastaRecord]:
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


def _parse_int(raw_value: str | None) -> int | None:
    if raw_value in {None, "", "."}:
        return None
    return int(float(raw_value))


def main() -> int:
    args = build_arg_parser().parse_args()
    write_te_fragments(
        args.sv_fasta,
        args.te_annotations,
        args.output,
        min_length=args.min_length,
        line_width=args.line_width,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())



