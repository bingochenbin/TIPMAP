from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

from tipmap.lib.fasta import iter_fasta

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "extract_te_fragments.py"
SPEC = importlib.util.spec_from_file_location("extract_te_fragments_script", SCRIPT_PATH)
assert SPEC is not None
extract_te_fragments_script = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = extract_te_fragments_script
assert SPEC.loader is not None
SPEC.loader.exec_module(extract_te_fragments_script)


class ExtractTeFragmentsScriptTests(unittest.TestCase):
    def test_extracts_te_fragments_from_annotation_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            fasta = base / "sv.fa"
            annotations = base / "te.tsv"
            output = base / "te.fa"
            fasta.write_text(
                ">sv1|alt|panpop|Chr1|1|1|INS|md5a\nAACCGGTTAACC\n",
                encoding="utf-8",
            )
            annotations.write_text(
                "seq_id\tmd5\tchrom\tstart\tend\tstrand\tsource\ttype\tfamily\tclass\tattributes\n"
                "sv1|alt|panpop|Chr1|1|1|INS|md5a\tmd5a\tChr1\t3\t8\t+\tEDTA\trepeat_region\tGypsy\tLTR\tID=x\n",
                encoding="utf-8",
            )

            count = extract_te_fragments_script.write_te_fragments(fasta, annotations, output)
            records = list(iter_fasta(output))

        self.assertEqual(count, 1)
        self.assertEqual(records[0].sequence, "CCGGTT")
        self.assertTrue(records[0].name.startswith("sv1|alt|panpop|Chr1|1|1|INS|md5a::te:3-8:"))
        self.assertIn("family=Gypsy", records[0].description or "")

    def test_reverse_complements_negative_strand_fragment(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            fasta = base / "sv.fa"
            annotations = base / "te.tsv"
            output = base / "te.fa"
            fasta.write_text(">seq1\nAACCGG\n", encoding="utf-8")
            annotations.write_text(
                "seq_id\tmd5\tchrom\tstart\tend\tstrand\tsource\ttype\tfamily\tclass\tattributes\n"
                "seq1\t\tChr1\t2\t5\t-\tEDTA\trepeat_region\tTIR\tDNA\tID=x\n",
                encoding="utf-8",
            )

            extract_te_fragments_script.write_te_fragments(fasta, annotations, output)
            record = next(iter_fasta(output))

        self.assertEqual(record.sequence, "CGGT")


if __name__ == "__main__":
    unittest.main()


