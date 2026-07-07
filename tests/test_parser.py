from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from tipmap.lib.parser import (
    detect_alt_kind,
    parse_info_field,
    parse_pairwise_vcf,
    parse_panpop_vcf,
    parse_vcf,
    parse_vcf_lines,
    parse_vcf_record_line,
)


class VCFParserTests(unittest.TestCase):
    def test_parse_info_field_supports_values_and_flags(self) -> None:
        info = parse_info_field("SVTYPE=INS;END=101256;SOMATIC;SVLEN=1234")
        self.assertEqual(info["SVTYPE"], "INS")
        self.assertEqual(info["END"], "101256")
        self.assertEqual(info["SVLEN"], "1234")
        self.assertEqual(info["SOMATIC"], "True")

    def test_detect_alt_kind_separates_real_and_symbolic_alt(self) -> None:
        self.assertEqual(detect_alt_kind("ACTG"), "sequence")
        self.assertEqual(detect_alt_kind("<INS>"), "symbolic")
        self.assertEqual(detect_alt_kind("N]Chr2:20]"), "breakend")
        self.assertEqual(detect_alt_kind("."), "missing")

    def test_parses_info_driven_insertion(self) -> None:
        records = parse_vcf_record_line(
            "Chr1A\t100023\tins-1\tN\tACTG\t60\tPASS\tSVTYPE=INS;END=101256;SVLEN=1234",
            genome="JM47",
        )

        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record.genome, "JM47")
        self.assertEqual(record.chrom, "Chr1A")
        self.assertEqual(record.pos, 100023)
        self.assertEqual(record.end, 101256)
        self.assertEqual(record.svtype, "INS")
        self.assertEqual(record.svlen, 1234)
        self.assertEqual(record.ref, "N")
        self.assertEqual(record.alt, "ACTG")
        self.assertTrue(record.has_real_alt_sequence)
        self.assertEqual(record.alt_kind, "sequence")
        self.assertEqual(record.source_id, "ins-1")
        self.assertEqual(record.qual, 60.0)
        self.assertEqual(record.filter, "PASS")

    def test_infers_svtype_and_svlen_from_sequence_alt(self) -> None:
        records = parse_vcf_record_line(
            "Chr2\t42\t.\tA\tATGCA\t.\t.\t.",
            genome="GenomeA",
        )

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].svtype, "INS")
        self.assertEqual(records[0].svlen, 4)
        self.assertEqual(records[0].end, 42)
        self.assertIsNone(records[0].source_id)
        self.assertIsNone(records[0].qual)
        self.assertIsNone(records[0].filter)

    def test_parses_symbolic_deletion_from_alt_and_end(self) -> None:
        records = parse_vcf_record_line(
            "Chr3\t10\tdel-1\tN\t<DEL>\t30\tPASS\tEND=25",
            genome="GenomeA",
        )

        self.assertEqual(records[0].svtype, "DEL")
        self.assertEqual(records[0].end, 25)
        self.assertEqual(records[0].svlen, -16)
        self.assertFalse(records[0].has_real_alt_sequence)
        self.assertTrue(records[0].is_symbolic_alt)

    def test_expands_multiallelic_ins_del_records_with_per_allele_info(self) -> None:
        records = parse_vcf_record_line(
            "Chr4\t100\tmulti-1\tN\t<DEL>,<INS>\t99\tPASS\tSVTYPE=DEL,INS;SVLEN=-50,80",
            genome="GenomeA",
        )

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].svtype, "DEL")
        self.assertEqual(records[0].svlen, -50)
        self.assertEqual(records[0].end, 149)
        self.assertEqual(records[0].allele_index, 0)
        self.assertEqual(records[1].svtype, "INS")
        self.assertEqual(records[1].svlen, 80)
        self.assertEqual(records[1].end, 100)
        self.assertEqual(records[1].allele_index, 1)

    def test_unsupported_svtypes_are_normalized_to_unknown(self) -> None:
        records = parse_vcf_record_line(
            "Chr4\t100\tdup-1\tN\t<DUP:TANDEM>\t99\tPASS\tSVTYPE=DUP;SVLEN=80",
            genome="GenomeA",
        )

        self.assertEqual(records[0].svtype, "UNK")
        self.assertEqual(records[0].svlen, 80)

    def test_parse_vcf_lines_uses_first_sample_as_genome(self) -> None:
        lines = [
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tJM47\n",
            "Chr1A\t100023\tins-1\tN\t<INS>\t60\tPASS\tSVTYPE=INS;END=101256;SVLEN=1234\tGT\t1/1\n",
        ]

        records = list(parse_vcf_lines(lines))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].genome, "JM47")
        self.assertEqual(records[0].sample_genotypes, {})

    def test_panpop_vcf_uses_file_stem_and_preserves_sample_genotypes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "panpop_SV.vcf"
            path.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCM104\tZ8425B",
                        "Chr1A\t703277\tChr1A_703277\tA\tAGCG\t.\tPASS\t.\tGT\t0/0\t1/1",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            record = list(parse_panpop_vcf(path))[0]

        self.assertEqual(record.genome, "panpop_SV")
        self.assertEqual(record.svtype, "INS")
        self.assertEqual(record.svlen, 3)
        self.assertEqual(record.sample_genotypes, {"CM104": "0/0", "Z8425B": "1/1"})
        self.assertFalse(record.carries_alt("CM104"))
        self.assertTrue(record.carries_alt("Z8425B"))

    def test_panpop_multiallelic_genotypes_are_allele_specific(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "panpop_multiallelic.vcf"
            path.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGenomeA\tGenomeB\tGenomeC",
                        "Chr1A\t1157425\tChr1A_1157425\tTCTTT\tA,ATCTC\t.\t.\tZLEN=5\tGT\t1/1\t2/2\t0/0",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            records = list(parse_panpop_vcf(path))

        self.assertEqual(len(records), 2)
        first_alt, second_alt = records
        self.assertEqual(first_alt.allele_index, 0)
        self.assertEqual(second_alt.allele_index, 1)
        self.assertTrue(first_alt.carries_alt("GenomeA"))
        self.assertFalse(first_alt.carries_alt("GenomeB"))
        self.assertFalse(first_alt.carries_alt("GenomeC"))
        self.assertFalse(second_alt.carries_alt("GenomeA"))
        self.assertTrue(second_alt.carries_alt("GenomeB"))
        self.assertFalse(second_alt.carries_alt("GenomeC"))

    def test_pairwise_vcf_preserves_assembly_coordinates_as_info(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "Z8425B.vcf"
            path.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tZ8425B",
                        (
                            "Chr1A\t1379602\t.\tA\tATTAA\t.\t.\t"
                            "ASM_Chr=1A;ASM_End=1399790;ASM_Start=1399357;ASM_Strand=+"
                            "\tGT:AD:DP:PL\t1:0,30,0:30:90,90"
                        ),
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            record = list(parse_pairwise_vcf(path))[0]

        self.assertEqual(record.genome, "Z8425B")
        self.assertEqual(record.svtype, "INS")
        self.assertEqual(record.svlen, 4)
        self.assertEqual(record.sample_genotypes, {})
        self.assertEqual(record.info["ASM_Start"], "1399357")
        self.assertEqual(record.info["ASM_End"], "1399790")

    def test_panpop_and_pairwise_wrappers_mark_vcf_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            panpop_path = Path(tmpdir) / "panpop.vcf"
            panpop_path.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                        "Chr1A\t100023\tpanpop-ins\tN\t<INS>\t60\tPASS\tSVTYPE=INS;SVLEN=1234",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            pairwise_path = Path(tmpdir) / "JM47.vcf"
            pairwise_path.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                        "Chr1A\t100023\tpairwise-ins\tN\tACTG\t60\tPASS\tSVTYPE=INS;SVLEN=3",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            panpop_record = list(parse_panpop_vcf(panpop_path, genome="JM47"))[0]
            pairwise_record = list(parse_pairwise_vcf(pairwise_path, genome="JM47"))[0]

        self.assertEqual(panpop_record.vcf_source, "panpop")
        self.assertEqual(panpop_record.alt_kind, "symbolic")
        self.assertFalse(panpop_record.has_real_alt_sequence)
        self.assertEqual(pairwise_record.vcf_source, "pairwise")
        self.assertEqual(pairwise_record.alt_kind, "sequence")
        self.assertTrue(pairwise_record.has_real_alt_sequence)

    def test_parse_vcf_reads_gzip_and_uses_filename_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "GenomeB.vcf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                handle.write("Chr1\t5\tdel-1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20\n")

            records = list(parse_vcf(path))

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].genome, "GenomeB")
        self.assertEqual(records[0].svtype, "DEL")
        self.assertEqual(records[0].svlen, -16)

    def test_parse_vcf_prefers_header_sample_over_filename(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "FileName.vcf"
            path.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tJM47",
                        "Chr1A\t100023\tins-1\tN\t<INS>\t60\tPASS\tSVTYPE=INS;SVLEN=1234\tGT\t1/1",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            records = list(parse_vcf(path))

        self.assertEqual(records[0].genome, "JM47")


if __name__ == "__main__":
    unittest.main()
