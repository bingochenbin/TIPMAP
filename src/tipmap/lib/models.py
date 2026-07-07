"""Typed data models shared by TIPMap modules."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Literal, Optional

VCFSource = Literal["panpop", "pairwise", "unknown"]
AltKind = Literal["sequence", "symbolic", "breakend", "missing"]
SequenceRole = Literal["ref", "alt"]


@dataclass(frozen=True)
class SVRecord:
    """A structural variant parsed from a VCF record."""

    genome: str
    chrom: str
    pos: int
    end: int
    svtype: str
    svlen: int
    ref: str
    alt: str
    source_id: Optional[str] = None
    qual: Optional[float] = None
    filter: Optional[str] = None
    info: Dict[str, str] = field(default_factory=dict)
    sample_genotypes: Dict[str, str] = field(default_factory=dict)
    allele_index: int = 0
    vcf_source: VCFSource = "unknown"
    alt_kind: AltKind = "sequence"

    @property
    def has_real_alt_sequence(self) -> bool:
        """Return ``True`` when ``alt`` contains an actual nucleotide sequence."""

        return self.alt_kind == "sequence"

    @property
    def is_symbolic_alt(self) -> bool:
        """Return ``True`` for symbolic VCF alleles such as ``<INS>``."""

        return self.alt_kind == "symbolic"

    def carries_alt(self, sample: str) -> bool:
        """Return whether ``sample`` carries this record's ALT allele."""

        genotype = self.sample_genotypes.get(sample)
        if genotype is None:
            return False
        alt_number = str(self.allele_index + 1)
        alleles = genotype.replace("|", "/").split("/")
        return any(allele == alt_number for allele in alleles)


@dataclass(frozen=True)
class ExtractedSequence:
    """A REF or ALT sequence extracted for one PanPop SV."""

    panpop_id: str
    role: SequenceRole
    genome: str
    chrom: str
    pos: int
    end: int
    svtype: str
    sequence: str
    md5: str
    source_id: Optional[str] = None
    allele_index: int = 0

    @property
    def length(self) -> int:
        """Return sequence length."""

        return len(self.sequence)
