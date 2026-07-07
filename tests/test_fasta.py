from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tipmap.lib.fasta import (
    FastaRecord,
    IndexedFasta,
    extract_interval,
    extract_interval_from_source,
    iter_fasta,
    read_fasta,
    reverse_complement,
    write_fasta,
)


class FastaTests(unittest.TestCase):
    def test_read_extract_and_write_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "ref.fa"
            output_path = Path(tmpdir) / "out.fa"
            input_path.write_text(
                ">Chr1 description\nAACCGGTTAACC\n>Chr2\nTTTT\n",
                encoding="utf-8",
            )

            sequences = read_fasta(input_path)
            self.assertEqual(sequences["Chr1"], "AACCGGTTAACC")
            self.assertEqual(extract_interval(sequences, "Chr1", 3, 6), "CCGG")
            self.assertEqual(extract_interval(sequences, "Chr1", 3, 6, flank=2), "AACCGGTT")

            write_fasta([FastaRecord(name="seq1", sequence="AACCGG", description="x=1")], output_path)
            records = list(iter_fasta(output_path))

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].name, "seq1")
        self.assertEqual(records[0].description, "x=1")
        self.assertEqual(records[0].sequence, "AACCGG")

    def test_extract_interval_rejects_invalid_coordinates(self) -> None:
        sequences = {"Chr1": "AACCGG"}
        with self.assertRaises(ValueError):
            extract_interval(sequences, "Chr1", 0, 2)
        with self.assertRaises(ValueError):
            extract_interval(sequences, "Chr1", 5, 2)
        with self.assertRaises(KeyError):
            extract_interval(sequences, "ChrX", 1, 2)


    def test_indexed_fasta_extracts_intervals_and_flanks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "ref.fa"
            input_path.write_text(
                ">Chr1 description\nAACCGG\nTTAACC\n>Chr2\nTTTT\n",
                encoding="utf-8",
            )

            indexed = IndexedFasta.open(input_path)

            self.assertIn("Chr1", indexed)
            self.assertEqual(indexed.extract_interval("Chr1", 3, 8), "CCGGTT")
            self.assertEqual(indexed.extract_interval("Chr1", 3, 8, flank=2), "AACCGGTTAA")
            self.assertEqual(extract_interval_from_source(indexed, "Chr2", 1, 4), "TTTT")


    def test_indexed_fasta_finds_stem_fai_next_to_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "ref.fa"
            primary_index = Path(str(input_path) + ".fai")
            stem_index = input_path.with_suffix(".fai")
            input_path.write_text(">Chr1\nAACCGG\n", encoding="utf-8")
            IndexedFasta.open(input_path)
            primary_index.replace(stem_index)

            indexed = IndexedFasta.open(input_path, index_mode="require")

            self.assertEqual(indexed.extract_interval("Chr1", 2, 5), "ACCG")
            self.assertFalse(primary_index.exists())
            self.assertTrue(stem_index.exists())

    def test_indexed_fasta_auto_generates_primary_fai_when_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "ref.fa"
            primary_index = Path(str(input_path) + ".fai")
            input_path.write_text(">Chr1\nAACCGG\n", encoding="utf-8")

            IndexedFasta.open(input_path, index_mode="auto")

            self.assertTrue(primary_index.exists())
    def test_indexed_fasta_require_mode_rejects_missing_index(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "ref.fa"
            input_path.write_text(">Chr1\nAACCGG\n", encoding="utf-8")

            with self.assertRaises(FileNotFoundError):
                IndexedFasta.open(input_path, index_mode="require")
    def test_reverse_complement(self) -> None:
        self.assertEqual(reverse_complement("AaCcGgTtNn"), "nNaAcCgGtT")


if __name__ == "__main__":
    unittest.main()



