from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "parse_edta_gff3.py"
SPEC = importlib.util.spec_from_file_location("parse_edta_gff3_script", SCRIPT_PATH)
assert SPEC is not None
parse_edta_gff3_script = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = parse_edta_gff3_script
assert SPEC.loader is not None
SPEC.loader.exec_module(parse_edta_gff3_script)


class ParseEdtaGff3ScriptTests(unittest.TestCase):
    def test_single_gff3_writes_standard_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            gff3 = base / "sample.EDTA.TEanno.gff3"
            output = base / "te.tsv"
            gff3.write_text(
                "##gff-version 3\n"
                "sv1|alt|GenomeA|Chr1A|1|1|INS|md5a\tEDTA\trepeat_region\t2\t20\t.\t+\t.\tID=x;Classification=LTR/Gypsy\n"
                "not\tenough\tfields\n",
                encoding="utf-8",
            )

            count = parse_edta_gff3_script.run_parse_workflow(gff3_files=[gff3], output=output)
            lines = output.read_text(encoding="utf-8").strip().splitlines()

        self.assertEqual(count, 1)
        self.assertEqual(
            lines[0],
            "seq_id\tmd5\tchrom\tstart\tend\tstrand\tsource\ttype\tfamily\tclass\tattributes",
        )
        self.assertIn("sv1|alt|GenomeA|Chr1A|1|1|INS|md5a\tmd5a\tChr1A\t2\t20", lines[1])
        self.assertIn("Gypsy\tLTR", lines[1])

    def test_multiple_gff3_inputs_are_merged_in_sorted_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            first = base / "b.gff3"
            second = base / "a.gff3"
            output = base / "te.tsv"
            first.write_text(
                "sv2|alt|GenomeA|Chr2A|1|1|INS|md5b\tEDTA\trepeat_region\t1\t5\t.\t+\t.\tID=b;Classification=DNA/TIR\n",
                encoding="utf-8",
            )
            second.write_text(
                "sv1|alt|GenomeA|Chr1A|1|1|INS|md5a\tEDTA\trepeat_region\t1\t5\t.\t+\t.\tID=a;Classification=LTR/Copia\n",
                encoding="utf-8",
            )

            count = parse_edta_gff3_script.run_parse_workflow(gff3_files=[first, second], output=output)
            lines = output.read_text(encoding="utf-8").strip().splitlines()

        self.assertEqual(count, 2)
        self.assertTrue(lines[1].startswith("sv1|alt"))
        self.assertTrue(lines[2].startswith("sv2|alt"))

    def test_edta_dir_discovers_only_top_level_gff3(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            edta_dir = base / "edta"
            nested = edta_dir / "Chr1A"
            nested.mkdir(parents=True)
            nested_gff3 = nested / "nested.EDTA.TEanno.gff3"
            top_gff3 = edta_dir / "result.EDTA.TEanno.gff3"
            output = base / "te.tsv"
            nested_gff3.write_text(
                "sv_nested|alt|GenomeA|Chr1A|1|1|INS|md5n\tEDTA\trepeat_region\t1\t4\t.\t-\t.\tName=helitron|Helitron\n",
                encoding="utf-8",
            )
            top_gff3.write_text(
                "sv1|alt|GenomeA|Chr1A|1|1|INS|md5a\tEDTA\trepeat_region\t1\t4\t.\t-\t.\tName=helitron|Helitron\n",
                encoding="utf-8",
            )

            paths = parse_edta_gff3_script.collect_gff3_paths([], [edta_dir])
            count = parse_edta_gff3_script.run_parse_workflow(edta_dirs=[edta_dir], output=output)
            lines = output.read_text(encoding="utf-8").strip().splitlines()

        self.assertEqual(paths, [top_gff3])
        self.assertEqual(count, 1)
        self.assertTrue(lines[1].startswith("sv1|alt"))

    def test_edta_dir_uses_first_nonempty_priority_pattern(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            edta_dir = base / "edta"
            edta_dir.mkdir()
            preferred = edta_dir / "preferred.EDTA.gff3"
            fallback = edta_dir / "fallback.gff3"
            preferred.write_text(
                "sv1|alt|GenomeA|Chr1A|1|1|INS|md5a\tEDTA\trepeat_region\t1\t4\t.\t+\t.\tClassification=LTR/Gypsy\n",
                encoding="utf-8",
            )
            fallback.write_text(
                "sv2|alt|GenomeA|Chr2A|1|1|INS|md5b\tEDTA\trepeat_region\t1\t4\t.\t+\t.\tClassification=DNA/TIR\n",
                encoding="utf-8",
            )

            paths = parse_edta_gff3_script.collect_gff3_paths([], [edta_dir])

        self.assertEqual(paths, [preferred])

    def test_main_requires_input_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "te.tsv"
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                with self.assertRaises(SystemExit) as error:
                    parse_edta_gff3_script.main(["--output", str(output)])

        self.assertEqual(error.exception.code, 2)

    def test_empty_edta_dir_errors(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            output = base / "te.tsv"
            with self.assertRaises(FileNotFoundError):
                parse_edta_gff3_script.run_parse_workflow(edta_dirs=[base], output=output)


if __name__ == "__main__":
    unittest.main()






