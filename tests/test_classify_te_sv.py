from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path
from typing import Sequence

from classify_te_sv import (
    TEHit,
    TEMetadata,
    classify_allele,
    iter_panpop_te_alleles,
    iter_te_fragments,
    process_vcf_record_line,
    read_blast_hits,
    read_te_metadata,
    run_classification_workflow,
    rewrite_genotype,
    union_interval_length,
    write_tip_outputs,
)
from tipmap.lib.fasta import iter_fasta
from tipmap.lib.utils import sequence_md5


class ClassifyTeSvTests(unittest.TestCase):
    def test_union_interval_length_merges_overlaps(self) -> None:
        self.assertEqual(union_interval_length([(10, 80), (60, 130), (200, 260)]), 182)

    def test_classify_allele_uses_union_coverage_and_weighted_identity(self) -> None:
        sequence = "A" * 100
        evidence = classify_allele(
            sequence,
            [TEHit(1, 30, 90.0, "Gypsy", "LTR"), TEHit(20, 70, 80.0, "Gypsy", "LTR")],
            min_te_coverage=0.60,
            min_identity=80.0,
            min_te_covered_bp=40,
        )

        self.assertTrue(evidence.is_te)
        self.assertEqual(evidence.te_covered_bp, 70)
        self.assertAlmostEqual(evidence.coverage, 0.70)
        self.assertAlmostEqual(evidence.weighted_identity or 0.0, (30 * 90.0 + 51 * 80.0) / 81)

    def test_rewrite_genotype_maps_removed_alt_to_reference(self) -> None:
        self.assertEqual(rewrite_genotype("1/1", {1: 1, 3: 2}), "1/1")
        self.assertEqual(rewrite_genotype("2/2", {1: 1, 3: 2}), "0/0")
        self.assertEqual(rewrite_genotype("3/3", {1: 1, 3: 2}), "2/2")
        self.assertEqual(rewrite_genotype("1|3", {1: 1, 3: 2}), "1|2")

    def test_read_blast_hits_uses_query_coordinates_and_identity(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            blast = base / "blast.tsv"
            allele_md5 = sequence_md5("A" * 100)
            blast.write_text(
                "%s\tte1\t87.5\t60\t0\t0\t5\t64\t1\t60\t1e-10\t120\n" % allele_md5,
                encoding="utf-8",
            )

            hits = read_blast_hits(
                blast,
                {"te1": TEMetadata(seq_id="te1", family="Gypsy", superfamily="LTR", start="1", end="60", strand="+", feature_type="repeat_region", attributes="ID=x")},
            )

        self.assertEqual(len(hits[allele_md5]), 1)
        self.assertEqual(hits[allele_md5][0].start, 5)
        self.assertEqual(hits[allele_md5][0].end, 64)
        self.assertEqual(hits[allele_md5][0].identity, 87.5)
        self.assertEqual(hits[allele_md5][0].family, "Gypsy")

    def test_read_blast_hits_maps_extracted_te_fragment_subject_to_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            blast = base / "blast.tsv"
            allele_md5 = sequence_md5("A" * 100)
            subject = "source|alt|panpop|Chr1|1|1|INS|abc::te:1-50:def"
            blast.write_text(
                "%s\t%s\t91\t50\t0\t0\t1\t50\t1\t50\t1e-10\t120\n" % (allele_md5, subject),
                encoding="utf-8",
            )

            hits = read_blast_hits(
                blast,
                {
                    "source|alt|panpop|Chr1|1|1|INS|abc": TEMetadata(
                        seq_id="source|alt|panpop|Chr1|1|1|INS|abc",
                        family="TIR",
                        superfamily="DNA",
                        start="1",
                        end="50",
                        strand="+",
                        feature_type="repeat_region",
                        attributes="ID=y",
                    )
                },
            )

        self.assertEqual(hits[allele_md5][0].family, "TIR")
        self.assertEqual(hits[allele_md5][0].superfamily, "DNA")
        self.assertEqual(hits[allele_md5][0].metadata.attributes if hits[allele_md5][0].metadata else "", "ID=y")

    def test_read_blast_hits_prefers_exact_te_fragment_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            te = base / "te.tsv"
            blast = base / "blast.tsv"
            seq_id = "source|alt|panpop|Chr1|1|1|INS|abc"
            allele_md5 = sequence_md5("A" * 120)
            te.write_text(
                "seq_id\tmd5\tchrom\tstart\tend\tstrand\tsource\ttype\tfamily\tsuperfamily\tattributes\n"
                "%s\tabc\tChr1\t1\t40\t+\tEDTA\trepeat_region\tGypsy\tLTR\tID=first\n"
                "%s\tabc\tChr1\t60\t100\t+\tEDTA\trepeat_region\tTIR\tDNA\tID=second\n" % (seq_id, seq_id),
                encoding="utf-8",
            )
            blast.write_text(
                "%s\t%s::te:1-40:frag1\t91\t40\t0\t0\t1\t40\t1\t40\t1e-10\t120\n"
                "%s\t%s::te:60-100:frag2\t92\t41\t0\t0\t60\t100\t1\t41\t1e-10\t120\n" % (allele_md5, seq_id, allele_md5, seq_id),
                encoding="utf-8",
            )

            metadata = read_te_metadata(te)
            hits = read_blast_hits(blast, metadata)[allele_md5]

        self.assertEqual([hit.family for hit in hits], ["Gypsy", "TIR"])
        self.assertEqual([hit.metadata.attributes if hit.metadata else "" for hit in hits], ["ID=first", "ID=second"])

    def test_process_multiallelic_ins_retains_only_te_alt_sequences(self) -> None:
        ref = "A"
        alt1 = "A" * 100
        alt2 = "C" * 100
        alt3 = "G" * 100
        hits_by_md5 = {
            sequence_md5(alt1): [TEHit(1, 70, 90.0, "Gypsy", "LTR")],
            sequence_md5(alt3): [TEHit(1, 80, 85.0, "Copia", "LTR")],
        }
        line = (
            "Chr1\t10\tsv1\t%s\t%s,%s,%s\t.\tPASS\t.\tGT\t1/1\t2/2\t3/3\t1/3"
            % (ref, alt1, alt2, alt3)
        )

        processed = process_vcf_record_line(line, hits_by_md5)

        self.assertIsNotNone(processed.line)
        fields = processed.line.split("\t") if processed.line is not None else []
        self.assertEqual(fields[4], "%s,%s" % (alt1, alt3))
        self.assertEqual(fields[9:], ["1/1", "0/0", "2/2", "1/2"])
        self.assertIn("TIP_TE_ALTS=1,3", fields[7])
        self.assertIn("TIP_RETAINED_ALTS=1,3", fields[7])

    def test_integrated_workflow_generates_fastas_runs_blast_and_writes_tip_vcf(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            panpop = base / "panpop.vcf"
            sv_fasta = base / "sv.fa"
            te_tsv = base / "te.tsv"
            out_vcf = base / "tip.vcf"
            report = base / "report.tsv"
            workdir = base / "work"
            alt = "A" * 100
            te_seq_id = "te_source|alt|GenomeA|Chr1|1|1|INS|md5te"
            panpop.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n"
                "Chr1\t10\tins1\tA\t%s\t.\tPASS\t.\tGT\t1/1\n" % alt,
                encoding="utf-8",
            )
            sv_fasta.write_text(">%s\n%s\n" % (te_seq_id, alt), encoding="utf-8")
            te_tsv.write_text(
                "seq_id\tmd5\tchrom\tstart\tend\tstrand\tsource\ttype\tfamily\tsuperfamily\tattributes\n"
                "%s\tmd5te\tChr1\t1\t70\t+\tEDTA\trepeat_region\tGypsy\tLTR\tID=x\n" % te_seq_id,
                encoding="utf-8",
            )

            def runner(command: Sequence[str], cwd: Path | None) -> subprocess.CompletedProcess[str]:
                query = Path(command[command.index("-query") + 1])
                out = Path(command[command.index("-out") + 1])
                query_record = next(iter_fasta(query))
                subject = "%s::te:1-70:fragmentmd5" % te_seq_id
                out.write_text(
                    "%s\t%s\t90\t70\t0\t0\t1\t70\t1\t70\t1e-20\t200\n" % (query_record.name, subject),
                    encoding="utf-8",
                )
                return subprocess.CompletedProcess(command, 0, "", "")

            record_count, report_count = run_classification_workflow(
                panpop_vcf=panpop,
                sv_fasta=sv_fasta,
                te_annotations=te_tsv,
                output_vcf=out_vcf,
                allele_report=report,
                workdir=workdir,
                runner=runner,
            )
            data = [line for line in out_vcf.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
            allele_fasta_exists = (workdir / "panpop_alleles.fa").is_file()
            te_fasta_exists = (workdir / "te_fragments.fa").is_file()
            blast_tsv_exists = (workdir / "panpop_allele_vs_te.tsv").is_file()
            report_text = report.read_text(encoding="utf-8")

        self.assertEqual(record_count, 1)
        self.assertEqual(report_count, 1)
        self.assertTrue(allele_fasta_exists)
        self.assertTrue(te_fasta_exists)
        self.assertTrue(blast_tsv_exists)
        self.assertEqual(data[0].split("\t")[9], "1/1")
        self.assertIn("TIP_TE_ALTS=1", data[0])
        self.assertIn("supporting_te_annotations", report_text)
        self.assertIn("seq_id=%s" % te_seq_id, report_text)
        self.assertIn("family=Gypsy", report_text)
        self.assertIn("attributes=ID=x", report_text)

    def test_write_tip_outputs_keeps_del_when_ref_is_te(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            ref = "T" * 100
            alt = "T"
            vcf = base / "panpop.vcf"
            blast = base / "blast.tsv"
            te = base / "te.tsv"
            out_vcf = base / "tip.vcf"
            report = base / "alleles.tsv"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
                "Chr1\t20\tdel1\t%s\t%s\t.\tPASS\t.\tGT\t0/0\t1/1\n" % (ref, alt),
                encoding="utf-8",
            )
            blast.write_text(
                "%s\tte1\t90\t70\t0\t0\t1\t70\t1\t70\t1e-20\t200\n" % sequence_md5(ref),
                encoding="utf-8",
            )
            te.write_text(
                "seq_id\tmd5\tchrom\tstart\tend\tstrand\tsource\ttype\tfamily\tsuperfamily\tattributes\n"
                "te1\t\tChr1\t1\t70\t+\tEDTA\trepeat_region\tTIR\tDNA\tID=x\n",
                encoding="utf-8",
            )

            record_count, report_count = write_tip_outputs(
                panpop_vcf=vcf,
                blast_tsv=blast,
                te_annotations=te,
                output_vcf=out_vcf,
                allele_report=report,
            )
            lines = out_vcf.read_text(encoding="utf-8").strip().splitlines()
            data = [line for line in lines if not line.startswith("#")]

        self.assertEqual(record_count, 1)
        self.assertEqual(report_count, 1)
        self.assertEqual(data[0].split("\t")[9:], ["0/0", "1/1"])
        self.assertIn("TIP_TE_REF=1", data[0])


if __name__ == "__main__":
    unittest.main()










