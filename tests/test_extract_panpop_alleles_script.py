from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

from tipmap.lib.fasta import iter_fasta

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "extract_panpop_alleles.py"
SPEC = importlib.util.spec_from_file_location("extract_panpop_alleles_script", SCRIPT_PATH)
assert SPEC is not None
extract_panpop_alleles_script = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = extract_panpop_alleles_script
assert SPEC.loader is not None
SPEC.loader.exec_module(extract_panpop_alleles_script)


class ExtractPanpopAllelesScriptTests(unittest.TestCase):
    def test_extracts_ins_alts_and_del_ref_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            vcf = base / "panpop.vcf"
            output = base / "alleles.fa"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n"
                "Chr1\t10\tins1\tA\t" + "A" * 50 + "," + "C" * 60 + "\t.\tPASS\t.\tGT\t1/2\n"
                "Chr2\t20\tdel1\t" + "G" * 70 + "\tG," + "T" * 5 + "\t.\tPASS\t.\tGT\t1/2\n",
                encoding="utf-8",
            )

            count = extract_panpop_alleles_script.write_panpop_te_alleles(vcf, output)
            records = list(iter_fasta(output))

        self.assertEqual(count, 3)
        self.assertEqual([record.name.split("|")[1] for record in records], ["alt", "alt", "ref"])
        self.assertEqual([len(record.sequence) for record in records], [50, 60, 70])
        self.assertIn("original_allele=1", records[0].description or "")
        self.assertIn("original_allele=2", records[1].description or "")
        self.assertIn("original_allele=0", records[2].description or "")

    def test_min_length_filters_short_panpop_alleles(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            vcf = base / "panpop.vcf"
            output = base / "alleles.fa"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "Chr1\t10\tins1\tA\t" + "A" * 30 + "," + "C" * 50 + "\t.\tPASS\t.\n",
                encoding="utf-8",
            )

            count = extract_panpop_alleles_script.write_panpop_te_alleles(vcf, output, min_length=40)
            records = list(iter_fasta(output))

        self.assertEqual(count, 1)
        self.assertEqual(records[0].sequence, "C" * 50)


if __name__ == "__main__":
    unittest.main()
