"""FASTA parsing, interval extraction, random access, and writing helpers."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
import textwrap
from typing import Dict, Optional, TextIO, Union

from tipmap.lib.models import ExtractedSequence

FastaPath = Union[str, Path]


@dataclass(frozen=True)
class FastaRecord:
    """A FASTA record."""

    name: str
    sequence: str
    description: Optional[str] = None


@dataclass(frozen=True)
class FastaIndexEntry:
    """One FASTA index entry compatible with the samtools ``.fai`` format."""

    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int


class IndexedFasta:
    """Random-access FASTA reader backed by a ``.fai``-style index."""

    def __init__(self, path: FastaPath, entries: Mapping[str, FastaIndexEntry]) -> None:
        self.path = Path(path)
        self.entries = dict(entries)

    @classmethod
    def open(cls, path: FastaPath, *, index_mode: str = "auto") -> "IndexedFasta":
        """Open a FASTA file with an existing or newly built index."""

        fasta_path = Path(path)
        if index_mode not in {"auto", "require", "rebuild"}:
            msg = "index_mode must be 'auto', 'require', or 'rebuild'"
            raise ValueError(msg)

        index_paths = _candidate_fai_paths(fasta_path)
        if index_mode != "rebuild":
            for index_path in index_paths:
                if index_path.exists():
                    return cls(fasta_path, _read_fai(index_path))
        if index_mode == "require":
            msg = "FASTA index not found beside %s; checked: %s" % (
                fasta_path,
                ", ".join(str(path) for path in index_paths),
            )
            raise FileNotFoundError(msg)

        index_path = index_paths[0]
        entries = _build_fai_entries(fasta_path)
        try:
            _write_fai(index_path, entries)
        except OSError:
            pass
        return cls(fasta_path, entries)

    def __contains__(self, chrom: object) -> bool:
        return chrom in self.entries

    def extract_interval(self, chrom: str, pos: int, end: int, *, flank: int = 0) -> str:
        """Extract a 1-based inclusive interval without loading the whole FASTA."""

        if chrom not in self.entries:
            msg = "Chromosome not found in FASTA: %s" % chrom
            raise KeyError(msg)
        if pos < 1:
            msg = "FASTA interval positions are 1-based and must be >= 1"
            raise ValueError(msg)
        if end < pos:
            msg = "FASTA interval end must be greater than or equal to pos"
            raise ValueError(msg)
        if flank < 0:
            msg = "flank must be >= 0"
            raise ValueError(msg)

        entry = self.entries[chrom]
        start = max(pos - flank, 1)
        stop = min(end + flank, entry.length)
        if stop < start:
            return ""

        chunks: list[str] = []
        remaining = stop - start + 1
        zero_based = start - 1
        with self.path.open("rb") as handle:
            while remaining > 0:
                line_index = zero_based // entry.line_bases
                column_index = zero_based % entry.line_bases
                take = min(remaining, entry.line_bases - column_index)
                byte_offset = entry.offset + line_index * entry.line_width + column_index
                handle.seek(byte_offset)
                chunks.append(handle.read(take).decode("ascii"))
                zero_based += take
                remaining -= take
        return "".join(chunks)


def iter_fasta(path: FastaPath) -> Iterator[FastaRecord]:
    """Yield FASTA records from ``path``."""

    name: Optional[str] = None
    description: Optional[str] = None
    chunks: list[str] = []
    with Path(path).open("rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield FastaRecord(name=name, sequence="".join(chunks), description=description)
                header = line[1:].strip()
                name, description = _split_header(header)
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield FastaRecord(name=name, sequence="".join(chunks), description=description)


def read_fasta(path: FastaPath) -> Dict[str, str]:
    """Read a FASTA file into a ``name -> sequence`` dictionary."""

    return {record.name: record.sequence for record in iter_fasta(path)}


def extract_interval(
    sequences: Mapping[str, str],
    chrom: str,
    pos: int,
    end: int,
    *,
    flank: int = 0,
) -> str:
    """Extract a 1-based inclusive interval from a FASTA sequence mapping."""

    if chrom not in sequences:
        msg = "Chromosome not found in FASTA: %s" % chrom
        raise KeyError(msg)
    if pos < 1:
        msg = "FASTA interval positions are 1-based and must be >= 1"
        raise ValueError(msg)
    if end < pos:
        msg = "FASTA interval end must be greater than or equal to pos"
        raise ValueError(msg)
    if flank < 0:
        msg = "flank must be >= 0"
        raise ValueError(msg)

    sequence = sequences[chrom]
    start_index = max(pos - flank - 1, 0)
    end_index = min(end + flank, len(sequence))
    return sequence[start_index:end_index]


def extract_interval_from_source(
    source: Mapping[str, str] | IndexedFasta,
    chrom: str,
    pos: int,
    end: int,
    *,
    flank: int = 0,
) -> str:
    """Extract an interval from either an in-memory FASTA mapping or an index."""

    if isinstance(source, IndexedFasta):
        return source.extract_interval(chrom, pos, end, flank=flank)
    return extract_interval(source, chrom, pos, end, flank=flank)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a nucleotide sequence."""

    complement = str.maketrans("ACGTRYMKSWBDHVNacgtrymkswbdhvn", "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn")
    return sequence.translate(complement)[::-1]


def write_fasta(
    records: Iterable[FastaRecord],
    path: FastaPath,
    *,
    line_width: int = 80,
) -> None:
    """Write FASTA records to ``path``."""

    with Path(path).open("wt", encoding="utf-8", newline="\n") as handle:
        write_fasta_records(records, handle, line_width=line_width)


def write_fasta_records(
    records: Iterable[FastaRecord],
    handle: TextIO,
    *,
    line_width: int = 80,
) -> None:
    """Write FASTA records to an open text handle."""

    for record in records:
        write_one_fasta_record(record, handle, line_width=line_width)


def write_one_fasta_record(
    record: FastaRecord,
    handle: TextIO,
    *,
    line_width: int = 80,
) -> None:
    """Write one FASTA record to an open text handle."""

    header = record.name
    if record.description:
        header += " " + record.description
    handle.write(">%s\n" % header)
    for line in textwrap.wrap(record.sequence, width=line_width) or [""]:
        handle.write("%s\n" % line)


def extracted_sequence_to_fasta_record(record: ExtractedSequence) -> FastaRecord:
    """Convert an extracted SV sequence to a FASTA record."""

    name = "%s|%s|%s|%s|%d|%d|%s|%s" % (
        record.panpop_id,
        record.role,
        record.genome,
        record.chrom,
        record.pos,
        record.end,
        record.svtype,
        record.md5,
    )
    description = "source_id=%s allele_index=%d length=%d" % (
        record.source_id or ".",
        record.allele_index,
        record.length,
    )
    return FastaRecord(name=name, sequence=record.sequence, description=description)


def _fai_path(path: Path) -> Path:
    return Path(str(path) + ".fai")


def _candidate_fai_paths(path: Path) -> tuple[Path, ...]:
    primary = _fai_path(path)
    secondary = path.with_suffix(".fai")
    if secondary == primary:
        return (primary,)
    return (primary, secondary)


def _read_fai(path: Path) -> Dict[str, FastaIndexEntry]:
    entries: Dict[str, FastaIndexEntry] = {}
    with path.open("rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                msg = "Invalid FASTA index line: %s" % line
                raise ValueError(msg)
            name = fields[0]
            entries[name] = FastaIndexEntry(
                name=name,
                length=int(fields[1]),
                offset=int(fields[2]),
                line_bases=int(fields[3]),
                line_width=int(fields[4]),
            )
    return entries


def _write_fai(path: Path, entries: Mapping[str, FastaIndexEntry]) -> None:
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        for entry in entries.values():
            handle.write(
                "%s\t%d\t%d\t%d\t%d\n"
                % (entry.name, entry.length, entry.offset, entry.line_bases, entry.line_width)
            )


def _build_fai_entries(path: Path) -> Dict[str, FastaIndexEntry]:
    entries: Dict[str, FastaIndexEntry] = {}
    current_name: Optional[str] = None
    current_length = 0
    sequence_offset = 0
    line_bases = 0
    line_width = 0

    def finish_record() -> None:
        if current_name is None:
            return
        if line_bases == 0 or line_width == 0:
            msg = "FASTA record has no sequence: %s" % current_name
            raise ValueError(msg)
        entries[current_name] = FastaIndexEntry(
            name=current_name,
            length=current_length,
            offset=sequence_offset,
            line_bases=line_bases,
            line_width=line_width,
        )

    with path.open("rb") as handle:
        while True:
            line_start = handle.tell()
            raw_line = handle.readline()
            if not raw_line:
                break
            if line_start == 0 and raw_line.startswith(b"\xef\xbb\xbf"):
                raw_line = raw_line[3:]
                line_start += 3
            if raw_line.startswith(b">"):
                finish_record()
                header = raw_line[1:].strip().decode("utf-8")
                current_name, _description = _split_header(header)
                current_length = 0
                sequence_offset = handle.tell()
                line_bases = 0
                line_width = 0
                continue
            if current_name is None:
                msg = "FASTA sequence encountered before first header"
                raise ValueError(msg)
            stripped = raw_line.rstrip(b"\r\n")
            if not stripped:
                continue
            bases = len(stripped)
            width = len(raw_line)
            if line_bases == 0:
                line_bases = bases
                line_width = width
            current_length += bases
        finish_record()
    return entries


def _split_header(header: str) -> tuple[str, Optional[str]]:
    if not header:
        msg = "FASTA header is empty"
        raise ValueError(msg)
    parts = header.split(maxsplit=1)
    if len(parts) == 1:
        return parts[0], None
    return parts[0], parts[1]


__all__ = [
    "FastaRecord",
    "FastaIndexEntry",
    "IndexedFasta",
    "iter_fasta",
    "read_fasta",
    "extract_interval",
    "extract_interval_from_source",
    "reverse_complement",
    "write_fasta",
    "write_fasta_records",
    "write_one_fasta_record",
    "extracted_sequence_to_fasta_record",
]


