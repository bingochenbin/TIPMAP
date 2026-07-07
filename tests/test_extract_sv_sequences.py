from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from extract_sv_sequences import (
    collect_extracted_sequences,
    collect_pairwise_records,
    collect_query_sequences_by_genome,
    read_pairwise_vcf_list,
)
from tipmap.lib.matcher import MatchCriteria
from tipmap.lib.models import SVRecord


class ExtractAltCliTests(unittest.TestCase):
    def test_read_pairwise_vcf_list_accepts_header_comments_and_relative_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            manifest = base / "pairwise_vcfs.tsv"
            manifest.write_text(
                "\n".join(
                    [
                        "# genome-to-reference VCF inputs",
                        "genome\tvcf_path",
                        "GenomeA\tGenomeA.vcf",
                        "GenomeB\tvcfs/GenomeB.vcf",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            records = read_pairwise_vcf_list(manifest)

        self.assertEqual([record.genome for record in records], ["GenomeA", "GenomeB"])
        self.assertEqual(records[0].path, base / "GenomeA.vcf")
        self.assertEqual(records[1].path, base / "vcfs" / "GenomeB.vcf")


    def test_read_pairwise_vcf_list_accepts_query_fasta_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            manifest = base / "pairwise_vcfs.tsv"
            manifest.write_text(
                "genome\tvcf_path\tquery_fasta_path\nGenomeA\tGenomeA.vcf\tGenomeA.fa\n",
                encoding="utf-8",
            )

            records = read_pairwise_vcf_list(manifest)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].query_fasta, base / "GenomeA.fa")

    def test_collect_query_sequences_by_genome_reads_manifest_fastas(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            query_fasta = base / "GenomeA.fa"
            query_fasta.write_text(">Chr1\nAACCGGTT\n", encoding="utf-8")
            manifest = base / "pairwise_vcfs.tsv"
            manifest.write_text("GenomeA\tGenomeA.vcf\tGenomeA.fa\n", encoding="utf-8")

            inputs = read_pairwise_vcf_list(manifest)
            sequences = collect_query_sequences_by_genome(inputs)

        self.assertEqual(sequences, {"GenomeA": {"Chr1": "AACCGGTT"}})
    def test_collect_pairwise_records_uses_manifest_genome_name(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            vcf_path = base / "not_the_genome_name.vcf"
            vcf_path.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                        "Chr1\t10\tins-1\tN\tNAACC\t.\tPASS\tSVTYPE=INS;SVLEN=4",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            manifest = base / "pairwise_vcfs.tsv"
            manifest.write_text("GenomeFromManifest\tnot_the_genome_name.vcf\n", encoding="utf-8")

            records = collect_pairwise_records([], pairwise_vcf_list=str(manifest))

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].genome, "GenomeFromManifest")
        self.assertEqual(records[0].alt, "NAACC")

    def test_collect_pairwise_records_parallel_matches_serial(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            first = base / "first.vcf"
            second = base / "second.vcf"
            first.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                        "Chr1\t10\tins-1\tN\tNAACC\t.\tPASS\tSVTYPE=INS;SVLEN=4",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            second.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                        "Chr1\t20\tins-2\tN\tNTTTT\t.\tPASS\tSVTYPE=INS;SVLEN=4",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            manifest = base / "pairwise_vcfs.tsv"
            manifest.write_text(
                "GenomeA\tfirst.vcf\nGenomeB\tsecond.vcf\n",
                encoding="utf-8",
            )

            serial = collect_pairwise_records([], pairwise_vcf_list=str(manifest), workers=1)
            parallel = collect_pairwise_records([], pairwise_vcf_list=str(manifest), workers=2)

        self.assertEqual(serial, parallel)
        self.assertEqual([record.genome for record in parallel], ["GenomeA", "GenomeB"])

    def test_collect_pairwise_records_rejects_invalid_worker_count(self) -> None:
        with self.assertRaises(ValueError):
            collect_pairwise_records(["example.vcf"], workers=0)

    def test_collect_extracted_sequences_parallel_matches_serial(self) -> None:
        panpop_records = [
            SVRecord(
                genome="panpop",
                chrom="Chr1",
                pos=2,
                end=2,
                svtype="INS",
                svlen=3,
                ref="A",
                alt="ATTT",
                source_id="ins-1",
                vcf_source="panpop",
                alt_kind="sequence",
                sample_genotypes={"GenomeA": "1/1"},
            ),
            SVRecord(
                genome="panpop",
                chrom="Chr1",
                pos=5,
                end=7,
                svtype="DEL",
                svlen=-3,
                ref="CCC",
                alt="C",
                source_id="del-1",
                vcf_source="panpop",
                alt_kind="sequence",
                sample_genotypes={"GenomeB": "1/1"},
            ),
        ]
        pairwise_records = [
            SVRecord(
                genome="GenomeA",
                chrom="Chr1",
                pos=2,
                end=2,
                svtype="INS",
                svlen=3,
                ref="A",
                alt="AGGG",
                vcf_source="pairwise",
                alt_kind="sequence",
            ),
            SVRecord(
                genome="GenomeB",
                chrom="Chr1",
                pos=5,
                end=7,
                svtype="DEL",
                svlen=-3,
                ref="CCC",
                alt="C",
                vcf_source="pairwise",
                alt_kind="sequence",
            ),
        ]
        reference_sequences = {"Chr1": "AACCCGGTT"}
        criteria = MatchCriteria(max_distance=0, max_length_ratio_difference=0.0)

        serial = collect_extracted_sequences(
            panpop_records,
            pairwise_records,
            reference_sequences,
            criteria=criteria,
            workers=1,
        )
        parallel = collect_extracted_sequences(
            panpop_records,
            pairwise_records,
            reference_sequences,
            criteria=criteria,
            workers=2,
        )

        self.assertEqual(serial, parallel)
        self.assertEqual([record.role for record in parallel], ["ref", "alt", "ref", "alt"])


    def test_collect_extracted_sequences_uses_query_alt_flank(self) -> None:
        panpop_records = [
            SVRecord(
                genome="panpop",
                chrom="Chr1",
                pos=2,
                end=2,
                svtype="INS",
                svlen=3,
                ref="A",
                alt="ATTT",
                source_id="ins-1",
                vcf_source="panpop",
                alt_kind="sequence",
                sample_genotypes={"GenomeA": "1/1"},
            )
        ]
        pairwise_records = [
            SVRecord(
                genome="GenomeA",
                chrom="Chr1",
                pos=2,
                end=2,
                svtype="INS",
                svlen=3,
                ref="A",
                alt="AGGG",
                vcf_source="pairwise",
                alt_kind="sequence",
                info={
                    "ASM_Chr": "Chr1",
                    "ASM_Start": "4",
                    "ASM_End": "6",
                    "ASM_Strand": "+",
                },
            )
        ]

        extracted = collect_extracted_sequences(
            panpop_records,
            pairwise_records,
            {"Chr1": "AACCGGTT"},
            criteria=MatchCriteria(max_distance=0, max_length_ratio_difference=0.0),
            alt_flank=2,
            query_sequences_by_genome={"GenomeA": {"Chr1": "AACCGGTTAA"}},
            workers=1,
        )

        self.assertEqual([record.role for record in extracted], ["ref", "alt"])
        self.assertEqual(extracted[1].sequence, "ACCGGTT")
    def test_collect_extracted_sequences_rejects_invalid_worker_count(self) -> None:
        with self.assertRaises(ValueError):
            collect_extracted_sequences(
                [],
                [],
                {},
                criteria=MatchCriteria(),
                workers=0,
            )


if __name__ == "__main__":
    unittest.main()



