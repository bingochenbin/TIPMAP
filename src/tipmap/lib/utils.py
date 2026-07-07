"""Small reusable helpers for TIPMap."""

from __future__ import annotations

import hashlib


def sequence_md5(sequence: str) -> str:
    """Return the MD5 digest used as the canonical sequence identifier."""

    return hashlib.md5(sequence.encode("utf-8")).hexdigest()


def gc_fraction(sequence: str) -> float:
    """Return GC fraction for a nucleotide sequence as a value from 0.0 to 1.0."""

    if not sequence:
        return 0.0
    gc_count = sum(1 for base in sequence.upper() if base in {"G", "C"})
    return gc_count / len(sequence)
