"""Extract PanPop allele sequences used for TE-SV/TIP classification."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterator

from tipmap.lib.fasta import FastaRecord, write_fasta_records
from tipmap.lib.matcher import panpop_sv_key
from tipmap.lib.parser import parse_panpop_vcf
from tipmap.lib.utils import sequence_md5


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panpop-vcf", required=True, help="Input PanPop merged SV VCF.")
    parser.add_argument("--output", required=True, help="Output PanPop allele FASTA path.")
    parser.add_argument(
        "--min-length",
        type=int,
        default=1,
        help="Minimum extracted allele sequence length. Default: 1.",
    )
    parser.add_argument(
        "--line-width",
        type=int,
        default=80,
        help="Output FASTA line width. Default: 80.",
    )
    return parser


def iter_panpop_te_alleles(panpop_vcf: str | Path, *, min_length: int = 1) -> Iterator[FastaRecord]:
    """Yield PanPop allele sequences for TE classification.

    INS records contribute their ALT allele sequence. DEL records contribute the
    REF allele sequence once per original PanPop record.
    """

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
        yield FastaRecord(name=name, sequence=sequence, description=description)


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


def main() -> int:
    args = build_arg_parser().parse_args()
    write_panpop_te_alleles(
        args.panpop_vcf,
        args.output,
        min_length=args.min_length,
        line_width=args.line_width,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
