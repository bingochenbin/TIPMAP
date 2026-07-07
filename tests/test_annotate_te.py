from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path
from typing import Sequence

from annotate_te import (
    SplitFasta,
    chromosome_from_tipmap_header,
    md5_from_tipmap_header,
    parse_edta_gff3,
    run_annotation_workflow,
    run_edta_by_chromosome,
    run_seqkit_rmdup,
    split_fasta_by_chromosome,
)
from tipmap.lib.fasta import iter_fasta


class AnnotateTeTests(unittest.TestCase):
    def test_header_helpers_parse_tipmap_fields(self) -> None:
        header = "sv1|alt|GenomeA|Chr1A|10|10|INS|abc123"

        self.assertEqual(chromosome_from_tipmap_header(header), "Chr1A")
        self.assertEqual(md5_from_tipmap_header(header), "abc123")

    def test_split_fasta_by_chromosome_uses_tipmap_header_chrom(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            fasta = base / "input.fa"
            fasta.write_text(
                ">sv1|alt|GenomeA|Chr1A|1|1|INS|md5a\nAAAA\n"
                ">sv2|ref|Ref|Chr2A|5|9|DEL|md5b\nCCCC\n",
                encoding="utf-8",
            )

            splits = split_fasta_by_chromosome(fasta, base / "splits")

            self.assertEqual([split.chrom for split in splits], ["Chr1A", "Chr2A"])
            chr1_records = list(iter_fasta(base / "splits" / "Chr1A.fa"))
            self.assertEqual(chr1_records[0].name, "sv1|alt|GenomeA|Chr1A|1|1|INS|md5a")

    def test_run_seqkit_rmdup_builds_sequence_based_command(self) -> None:
        calls: list[tuple[list[str], Path | None]] = []

        def runner(command: Sequence[str], cwd: Path | None) -> subprocess.CompletedProcess[str]:
            calls.append((list(command), cwd))
            return subprocess.CompletedProcess(command, 0, "", "")

        run_seqkit_rmdup("in.fa", "out.fa", seqkit="seqkit-bin", runner=runner)

        self.assertEqual(calls, [(["seqkit-bin", "rmdup", "-s", "-o", "out.fa", "in.fa"], None)])

    def test_run_edta_by_chromosome_builds_commands(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            split_a = base / "Chr1.fa"
            split_b = base / "Chr2.fa"
            split_a.write_text(">a\nAAAA\n", encoding="utf-8")
            split_b.write_text(">b\nCCCC\n", encoding="utf-8")
            calls: list[list[str]] = []

            def runner(command: Sequence[str], cwd: Path | None) -> subprocess.CompletedProcess[str]:
                calls.append(list(command))
                assert cwd is not None
                gff = cwd / "result.EDTA.TEanno.gff3"
                gff.write_text("##gff-version 3\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, "", "")

            results = run_edta_by_chromosome(
                [SplitFasta("Chr1", split_a, 1), SplitFasta("Chr2", split_b, 1)],
                base / "edta",
                edta="EDTA.pl",
                species="others",
                edta_threads=2,
                chrom_workers=2,
                extra_args=["--sensitive 1", "--curatedlib curated.fa"],
                runner=runner,
            )

        self.assertEqual(len(calls), 2)
        self.assertIn("--genome", calls[0])
        self.assertIn("--threads", calls[0])
        self.assertIn("--sensitive", calls[0])
        self.assertIn("1", calls[0])
        self.assertIn("--curatedlib", calls[0])
        self.assertIn("curated.fa", calls[0])
        self.assertEqual(len(results), 2)
        self.assertTrue(all(result.gff3 is not None for result in results))

    def test_run_edta_ignores_duplicate_anno_one(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            split = base / "Chr1.fa"
            split.write_text(">a\nAAAA\n", encoding="utf-8")
            calls: list[list[str]] = []

            def runner(command: Sequence[str], cwd: Path | None) -> subprocess.CompletedProcess[str]:
                calls.append(list(command))
                assert cwd is not None
                (cwd / "result.EDTA.TEanno.gff3").write_text("##gff-version 3\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, "", "")

            run_edta_by_chromosome(
                [SplitFasta("Chr1", split, 1)],
                base / "edta",
                extra_args=["--anno 1", "--anno=1", "--curatedlib curated.fa"],
                runner=runner,
            )

        self.assertEqual(calls[0].count("--anno"), 1)
        self.assertIn("--curatedlib", calls[0])
        self.assertIn("curated.fa", calls[0])

    def test_run_edta_warns_and_ignores_anno_zero(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            split = base / "Chr1.fa"
            split.write_text(">a\nAAAA\n", encoding="utf-8")
            calls: list[list[str]] = []

            def runner(command: Sequence[str], cwd: Path | None) -> subprocess.CompletedProcess[str]:
                calls.append(list(command))
                assert cwd is not None
                (cwd / "result.EDTA.TEanno.gff3").write_text("##gff-version 3\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, "", "")

            with self.assertLogs("annotate_te", level="WARNING") as logs:
                run_edta_by_chromosome(
                    [SplitFasta("Chr1", split, 1)],
                    base / "edta",
                    extra_args=["--anno 0"],
                    runner=runner,
                )

        self.assertEqual(calls[0].count("--anno"), 1)
        anno_index = calls[0].index("--anno")
        self.assertEqual(calls[0][anno_index + 1], "1")
        self.assertTrue(any("--anno 0" in message for message in logs.output))

    def test_parse_edta_gff3_normalizes_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            gff = Path(tmpdir) / "anno.gff3"
            gff.write_text(
                "##gff-version 3\n"
                "sv1|alt|GenomeA|Chr1A|1|1|INS|md5a\tEDTA\trepeat_region\t2\t20\t.\t+\t.\tID=x;Classification=LTR/Gypsy\n",
                encoding="utf-8",
            )

            rows = list(parse_edta_gff3(gff))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].seq_id, "sv1|alt|GenomeA|Chr1A|1|1|INS|md5a")
        self.assertEqual(rows[0].md5, "md5a")
        self.assertEqual(rows[0].chrom, "Chr1A")
        self.assertEqual(rows[0].family, "Gypsy")
        self.assertEqual(rows[0].superfamily, "LTR")

    def test_run_annotation_workflow_with_fake_runner(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            fasta = base / "input.fa"
            output = base / "te.tsv"
            workdir = base / "work"
            fasta.write_text(
                ">sv1|alt|GenomeA|Chr1A|1|1|INS|md5a\nAAAA\n"
                ">sv2|alt|GenomeA|Chr2A|1|1|INS|md5b\nCCCC\n",
                encoding="utf-8",
            )

            def runner(command: Sequence[str], cwd: Path | None) -> subprocess.CompletedProcess[str]:
                if command[1:3] == ["rmdup", "-s"]:
                    out_path = Path(command[command.index("-o") + 1])
                    out_path.write_text(fasta.read_text(encoding="utf-8"), encoding="utf-8")
                else:
                    assert cwd is not None
                    genome = Path(command[command.index("--genome") + 1])
                    first_record = next(iter_fasta(genome))
                    (cwd / "result.EDTA.TEanno.gff3").write_text(
                        "%s\tEDTA\trepeat_region\t1\t4\t.\t+\t.\tID=x;Classification=DNA/TIR\n"
                        % first_record.name,
                        encoding="utf-8",
                    )
                return subprocess.CompletedProcess(command, 0, "", "")

            count = run_annotation_workflow(
                fasta,
                output,
                workdir=workdir,
                seqkit="seqkit",
                edta="EDTA.pl",
                chrom_workers=2,
                runner=runner,
                dry_run=False,
            )

            lines = output.read_text(encoding="utf-8").strip().splitlines()

        self.assertEqual(count, 2)
        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[0].startswith("seq_id\tmd5\tchrom"))
        self.assertIn("DNA", lines[1])


if __name__ == "__main__":
    unittest.main()


