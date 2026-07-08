"""Run seqkit de-duplication and EDTA TE annotation for TIPMap sequences."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
import logging
import re
import shlex
import shutil
import subprocess
from typing import Callable, Iterable, Iterator, Sequence, TextIO

from tipmap.lib.fasta import iter_fasta, write_fasta_records

LOGGER = logging.getLogger(__name__)
CommandRunner = Callable[[Sequence[str], Path | None], subprocess.CompletedProcess[str]]


@dataclass(frozen=True)
class SplitFasta:
    """A chromosome-specific FASTA file produced from the de-duplicated input."""

    chrom: str
    path: Path
    count: int


@dataclass(frozen=True)
class EDTARunResult:
    """Files produced by one EDTA run."""

    chrom: str
    fasta: Path
    workdir: Path
    gff3: Path | None


@dataclass(frozen=True)
class TEAnnotationRow:
    """One parsed TE annotation row for TSV export."""

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
    parser.add_argument("fasta", help="Input FASTA produced by extract_sv_sequences.py.")
    parser.add_argument("--output", required=True, help="Final TE annotation TSV path.")
    parser.add_argument(
        "--workdir",
        required=True,
        help="Working directory for de-duplicated FASTA, chromosome splits, and EDTA outputs.",
    )
    parser.add_argument("--seqkit", default="seqkit", help="seqkit executable path.")
    parser.add_argument("--edta", default="EDTA.pl", help="EDTA.pl executable path.")
    parser.add_argument("--species", default="others", help="EDTA --species value.")
    parser.add_argument("--edta-threads", type=int, default=4, help="Threads passed to each EDTA run.")
    parser.add_argument(
        "--chrom-workers",
        type=int,
        default=1,
        help="Number of chromosome-specific EDTA jobs to run in parallel.",
    )
    parser.add_argument(
        "--fasta-line-width",
        type=int,
        default=80,
        help="Line width for chromosome-specific FASTA files.",
    )
    parser.add_argument(
        "--edta-arg",
        action="append",
        default=[],
        help="Extra EDTA.pl argument string. Repeat for multiple argument groups.",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep working files. Currently work files are kept by default for reproducibility.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print planned commands without running them.")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return parser


def run_seqkit_rmdup(
    input_fasta: str | Path,
    output_fasta: str | Path,
    *,
    seqkit: str = "seqkit",
    runner: CommandRunner | None = None,
    dry_run: bool = False,
) -> Path:
    """Run ``seqkit rmdup`` by sequence and return the de-duplicated FASTA path."""

    output_path = Path(output_fasta)
    command = [seqkit, "rmdup", "-s", "-o", str(output_path), str(input_fasta)]
    _run_command(command, runner=runner, dry_run=dry_run)
    return output_path


def split_fasta_by_chromosome(
    fasta: str | Path,
    output_dir: str | Path,
    *,
    line_width: int = 80,
) -> list[SplitFasta]:
    """Split a TIPMap FASTA by chromosome parsed from the FASTA header."""

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    handles: dict[str, TextIO] = {}
    counts: dict[str, int] = {}
    paths: dict[str, Path] = {}
    try:
        for record in iter_fasta(fasta):
            chrom = chromosome_from_tipmap_header(record.name)
            if chrom not in handles:
                split_path = output_path / (sanitize_filename(chrom) + ".fa")
                paths[chrom] = split_path
                handles[chrom] = split_path.open("wt", encoding="utf-8", newline="\n")
                counts[chrom] = 0
            write_fasta_records([record], handles[chrom], line_width=line_width)
            counts[chrom] += 1
    finally:
        for handle in handles.values():
            handle.close()
    return [SplitFasta(chrom=chrom, path=paths[chrom], count=counts[chrom]) for chrom in sorted(paths)]


def chromosome_from_tipmap_header(header_name: str) -> str:
    """Extract chromosome from TIPMap FASTA record name."""

    fields = header_name.split("|")
    if len(fields) >= 4 and fields[3]:
        return fields[3]
    return "unknown"


def sanitize_filename(value: str) -> str:
    """Return a filesystem-safe filename component."""

    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return sanitized or "unknown"


def append_edta_extra_args(command: list[str], extra_args: Sequence[str]) -> None:
    """Append user-provided EDTA arguments while keeping annotation enabled."""

    for extra in extra_args:
        tokens = shlex.split(extra)
        index = 0
        while index < len(tokens):
            token = tokens[index]
            if token == "--anno":
                value = tokens[index + 1] if index + 1 < len(tokens) else ""
                _log_ignored_anno_value(value)
                index += 2 if index + 1 < len(tokens) else 1
                continue
            if token.startswith("--anno="):
                _log_ignored_anno_value(token.split("=", 1)[1])
                index += 1
                continue
            command.append(token)
            index += 1


def _log_ignored_anno_value(value: str) -> None:
    if value == "1":
        LOGGER.debug("Ignoring duplicate EDTA --anno 1 argument; TIPMap always uses --anno 1.")
    elif value == "0":
        LOGGER.warning("Ignoring EDTA --anno 0 argument; TIPMap requires --anno 1 for TE annotation export.")
    else:
        LOGGER.warning("Ignoring unsupported EDTA --anno value %r; TIPMap requires --anno 1.", value)


def run_edta_for_split(
    split: SplitFasta,
    workdir: str | Path,
    *,
    edta: str = "EDTA.pl",
    species: str = "others",
    edta_threads: int = 4,
    extra_args: Sequence[str] = (),
    runner: CommandRunner | None = None,
    dry_run: bool = False,
) -> EDTARunResult:
    """Run EDTA for one chromosome-specific FASTA split."""

    chrom_workdir = Path(workdir) / sanitize_filename(split.chrom)
    chrom_workdir.mkdir(parents=True, exist_ok=True)
    command = [
        edta,
        "--genome",
        str(split.path.resolve()),
        "--species",
        species,
        "--threads",
        str(edta_threads),
        "--anno",
        "1",
    ]
    append_edta_extra_args(command, extra_args)
    _run_command(command, cwd=chrom_workdir, runner=runner, dry_run=dry_run)
    return EDTARunResult(
        chrom=split.chrom,
        fasta=split.path,
        workdir=chrom_workdir,
        gff3=None if dry_run else find_edta_gff3(chrom_workdir, split.path),
    )


def run_edta_by_chromosome(
    splits: Sequence[SplitFasta],
    workdir: str | Path,
    *,
    edta: str = "EDTA.pl",
    species: str = "others",
    edta_threads: int = 4,
    chrom_workers: int = 1,
    extra_args: Sequence[str] = (),
    runner: CommandRunner | None = None,
    dry_run: bool = False,
) -> list[EDTARunResult]:
    """Run chromosome-specific EDTA jobs, optionally in parallel."""

    if chrom_workers < 1:
        msg = "chrom_workers must be >= 1"
        raise ValueError(msg)
    if edta_threads < 1:
        msg = "edta_threads must be >= 1"
        raise ValueError(msg)

    def submit(split: SplitFasta) -> EDTARunResult:
        return run_edta_for_split(
            split,
            workdir,
            edta=edta,
            species=species,
            edta_threads=edta_threads,
            extra_args=extra_args,
            runner=runner,
            dry_run=dry_run,
        )

    if chrom_workers == 1 or len(splits) <= 1:
        return [submit(split) for split in splits]
    with ThreadPoolExecutor(max_workers=min(chrom_workers, len(splits))) as executor:
        return list(executor.map(submit, splits))


def find_edta_gff3(workdir: str | Path, split_fasta: str | Path) -> Path | None:
    """Find the most likely top-level EDTA GFF3 output for one split."""

    split_path = Path(split_fasta)
    search_dirs = [Path(workdir)]
    if split_path.parent not in search_dirs:
        search_dirs.append(split_path.parent)

    for search_dir in search_dirs:
        found = _find_top_level_edta_gff3(search_dir, split_path.name)
        if found is not None:
            return found
    return None


def _find_top_level_edta_gff3(directory: Path, split_stem: str) -> Path | None:
    for pattern in ("*EDTA.TEanno.gff3", "*TEanno*.gff3", "*EDTA*.gff3", "*.gff3"):
        candidates = sorted(path for path in directory.glob(pattern) if path.is_file())
        if candidates:
            break
    else:
        return None
    preferred = [path for path in candidates if split_stem in path.name]
    return preferred[0] if preferred else candidates[0]


def parse_edta_gff3(path: str | Path) -> Iterator[TEAnnotationRow]:
    """Parse EDTA GFF3 into normalized TSV rows."""

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
    """Parse GFF3 attributes into a string dictionary."""

    attributes: dict[str, str] = {}
    for item in raw_attributes.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes


def classify_te_attributes(attributes: dict[str, str], fallback: str) -> tuple[str, str]:
    """Extract family and superfamily from common EDTA GFF3 attribute keys."""

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


def md5_from_tipmap_header(header_name: str) -> str:
    """Extract the MD5 field from TIPMap FASTA record names."""

    fields = header_name.split("|")
    if len(fields) >= 8:
        return fields[7]
    return ""


def write_annotation_tsv(rows: Iterable[TEAnnotationRow], output: str | Path) -> int:
    """Write TE annotation rows to a normalized TSV file."""

    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with output_path.open("wt", encoding="utf-8", newline="\n") as handle:
        handle.write("seq_id\tmd5\tchrom\tstart\tend\tstrand\tsource\ttype\tfamily\tsuperfamily\tattributes\n")
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


def run_annotation_workflow(
    fasta: str | Path,
    output: str | Path,
    *,
    workdir: str | Path | None = None,
    seqkit: str = "seqkit",
    edta: str = "EDTA.pl",
    species: str = "others",
    edta_threads: int = 4,
    chrom_workers: int = 1,
    fasta_line_width: int = 80,
    edta_args: Sequence[str] = (),
    runner: CommandRunner | None = None,
    dry_run: bool = False,
) -> int:
    """Run seqkit de-duplication, chromosome-split EDTA, and TSV export."""

    input_fasta = Path(fasta)
    output_path = Path(output)
    work_path = Path(workdir) if workdir is not None else output_path.with_suffix(".work")
    work_path.mkdir(parents=True, exist_ok=True)
    if not dry_run and runner is None:
        _require_tool(seqkit)
        _require_tool(edta)

    dedup_fasta = work_path / (input_fasta.stem + ".rmdup.fa")
    run_seqkit_rmdup(input_fasta, dedup_fasta, seqkit=seqkit, runner=runner, dry_run=dry_run)
    splits = split_fasta_by_chromosome(dedup_fasta, work_path / "by_chrom", line_width=fasta_line_width)
    LOGGER.info("Created %d chromosome FASTA splits", len(splits))
    results = run_edta_by_chromosome(
        splits,
        work_path / "edta",
        edta=edta,
        species=species,
        edta_threads=edta_threads,
        chrom_workers=chrom_workers,
        extra_args=edta_args,
        runner=runner,
        dry_run=dry_run,
    )
    if dry_run:
        return 0
    rows = _iter_annotation_rows(results)
    return write_annotation_tsv(rows, output_path)


def _iter_annotation_rows(results: Iterable[EDTARunResult]) -> Iterator[TEAnnotationRow]:
    for result in results:
        if result.gff3 is None:
            LOGGER.warning("No EDTA GFF3 found for chromosome %s in %s", result.chrom, result.workdir)
            continue
        yield from parse_edta_gff3(result.gff3)


def _run_command(
    command: Sequence[str],
    *,
    cwd: Path | None = None,
    runner: CommandRunner | None = None,
    dry_run: bool = False,
) -> subprocess.CompletedProcess[str]:
    LOGGER.info("Running: %s", " ".join(command))
    if dry_run:
        return subprocess.CompletedProcess(command, 0, "", "")
    active_runner = runner or _default_runner
    result = active_runner(command, cwd)
    if result.returncode != 0:
        msg = "Command failed with exit code %d: %s" % (result.returncode, " ".join(command))
        raise RuntimeError(msg)
    return result


def _default_runner(command: Sequence[str], cwd: Path | None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(command),
        cwd=None if cwd is None else str(cwd),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _require_tool(executable: str) -> None:
    if Path(executable).exists():
        return
    if shutil.which(executable) is None:
        msg = "Required executable not found on PATH: %s" % executable
        raise FileNotFoundError(msg)


def main() -> int:
    args = build_arg_parser().parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")
    count = run_annotation_workflow(
        args.fasta,
        args.output,
        workdir=args.workdir,
        seqkit=args.seqkit,
        edta=args.edta,
        species=args.species,
        edta_threads=args.edta_threads,
        chrom_workers=args.chrom_workers,
        fasta_line_width=args.fasta_line_width,
        edta_args=args.edta_arg,
        dry_run=args.dry_run,
    )
    if args.dry_run:
        LOGGER.info("Dry run completed.")
    else:
        LOGGER.info("Wrote %d TE annotation rows to %s", count, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())







