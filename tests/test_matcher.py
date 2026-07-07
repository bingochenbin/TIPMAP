from __future__ import annotations

import unittest

from tipmap.lib.matcher import (
    MatchCriteria,
    PairwiseRecordIndex,
    choose_representative_alt,
    extract_alt_sequence,
    find_matching_alt_records,
    iter_extracted_sequences,
)
from tipmap.lib.models import ExtractedSequence, SVRecord
from tipmap.lib.utils import sequence_md5


class MatcherTests(unittest.TestCase):
    def test_matches_pairwise_true_alt_to_panpop_symbolic_sv(self) -> None:
        panpop = SVRecord(
            genome="JM47",
            chrom="Chr1A",
            pos=100,
            end=110,
            svtype="INS",
            svlen=10,
            ref="N",
            alt="<INS>",
            source_id="panpop-1",
            vcf_source="panpop",
            alt_kind="symbolic",
        )
        pairwise = SVRecord(
            genome="JM47",
            chrom="Chr1A",
            pos=103,
            end=111,
            svtype="INS",
            svlen=9,
            ref="N",
            alt="AACCGGTTA",
            source_id="pairwise-1",
            vcf_source="pairwise",
            alt_kind="sequence",
        )

        matches = find_matching_alt_records(
            panpop,
            [pairwise],
            criteria=MatchCriteria(max_distance=5, max_length_ratio_difference=0.2),
        )

        self.assertEqual(matches, [pairwise])

    def test_iter_extracted_sequences_yields_ref_and_all_matched_alt(self) -> None:
        panpop = SVRecord(
            genome="JM47",
            chrom="Chr1A",
            pos=2,
            end=5,
            svtype="DEL",
            svlen=-4,
            ref="ACCG",
            alt="<DEL>",
            source_id="panpop-del",
            vcf_source="panpop",
            alt_kind="symbolic",
            sample_genotypes={"GenomeQ": "1/1", "GenomeR": "1/1"},
        )
        pairwise_1 = SVRecord(
            genome="GenomeQ",
            chrom="Chr1A",
            pos=2,
            end=5,
            svtype="DEL",
            svlen=-4,
            ref="ACCG",
            alt="A",
            source_id="pairwise-del",
            vcf_source="pairwise",
            alt_kind="sequence",
        )
        pairwise_2 = SVRecord(
            genome="GenomeR",
            chrom="Chr1A",
            pos=2,
            end=5,
            svtype="DEL",
            svlen=-4,
            ref="ACCG",
            alt="AT",
            source_id="pairwise-del-2",
            vcf_source="pairwise",
            alt_kind="sequence",
        )

        extracted = list(
            iter_extracted_sequences(
                [panpop],
                [pairwise_1, pairwise_2],
                {"Chr1A": "AACCGGTT"},
                criteria=MatchCriteria(max_distance=0, max_length_ratio_difference=1.0),
            )
        )

        self.assertEqual([record.role for record in extracted], ["ref", "alt", "alt"])
        self.assertEqual(extracted[0].sequence, "ACCG")
        self.assertEqual(extracted[0].panpop_id, "panpop-del")
        self.assertEqual(extracted[1].sequence, "A")
        self.assertEqual(extracted[1].md5, sequence_md5("A"))
        self.assertEqual(extracted[2].sequence, "AT")
        self.assertEqual(extracted[2].genome, "GenomeR")

    def test_panpop_reference_genotype_suppresses_pairwise_alt_extraction(self) -> None:
        panpop = SVRecord(
            genome="panpop",
            chrom="Chr1A",
            pos=10,
            end=10,
            svtype="INS",
            svlen=4,
            ref="A",
            alt="ATTTT",
            source_id="panpop-ins",
            vcf_source="panpop",
            alt_kind="sequence",
            sample_genotypes={"GenomeRefLike": "0/0", "GenomeCarrier": "1/1"},
        )
        ref_like_pairwise = SVRecord(
            genome="GenomeRefLike",
            chrom="Chr1A",
            pos=10,
            end=10,
            svtype="INS",
            svlen=4,
            ref="A",
            alt="ATTTT",
            source_id="pairwise-ref-like",
            vcf_source="pairwise",
            alt_kind="sequence",
        )
        carrier_pairwise = SVRecord(
            genome="GenomeCarrier",
            chrom="Chr1A",
            pos=10,
            end=10,
            svtype="INS",
            svlen=4,
            ref="A",
            alt="AGGGG",
            source_id="pairwise-carrier",
            vcf_source="pairwise",
            alt_kind="sequence",
        )

        matches = find_matching_alt_records(
            panpop,
            [ref_like_pairwise, carrier_pairwise],
            criteria=MatchCriteria(max_distance=0, max_length_ratio_difference=0.0),
        )

        self.assertEqual(matches, [carrier_pairwise])

    def test_pairwise_record_index_uses_genome_chrom_type_and_position_window(self) -> None:
        panpop = SVRecord(
            genome="panpop",
            chrom="Chr1A",
            pos=1000,
            end=1000,
            svtype="INS",
            svlen=5,
            ref="A",
            alt="AAAAAA",
            source_id="panpop-ins",
            vcf_source="panpop",
            alt_kind="sequence",
            sample_genotypes={"GenomeA": "1/1", "GenomeB": "0/0"},
        )
        matching = SVRecord(
            genome="GenomeA",
            chrom="Chr1A",
            pos=1003,
            end=1002,
            svtype="INS",
            svlen=5,
            ref="A",
            alt="ACCCCC",
            source_id="match",
            vcf_source="pairwise",
            alt_kind="sequence",
        )
        wrong_genome = SVRecord(
            genome="GenomeB",
            chrom="Chr1A",
            pos=1000,
            end=1000,
            svtype="INS",
            svlen=5,
            ref="A",
            alt="AGGGGG",
            source_id="wrong-genome",
            vcf_source="pairwise",
            alt_kind="sequence",
        )
        wrong_chrom = SVRecord(
            genome="GenomeA",
            chrom="Chr2A",
            pos=1000,
            end=1000,
            svtype="INS",
            svlen=5,
            ref="A",
            alt="ATTTTT",
            source_id="wrong-chrom",
            vcf_source="pairwise",
            alt_kind="sequence",
        )
        distant = SVRecord(
            genome="GenomeA",
            chrom="Chr1A",
            pos=2000,
            end=2000,
            svtype="INS",
            svlen=5,
            ref="A",
            alt="ATAAAA",
            source_id="distant",
            vcf_source="pairwise",
            alt_kind="sequence",
        )
        pairwise_records = [wrong_genome, wrong_chrom, matching, distant]
        index = PairwiseRecordIndex(pairwise_records)

        matches = find_matching_alt_records(
            panpop,
            index,
            criteria=MatchCriteria(max_distance=5, max_length_ratio_difference=0.0),
        )

        self.assertEqual(matches, [matching])

    def test_multiallelic_panpop_genotype_matches_only_its_alt_allele(self) -> None:
        first_alt = SVRecord(
            genome="panpop",
            chrom="Chr1A",
            pos=100,
            end=104,
            svtype="DEL",
            svlen=-5,
            ref="AAAAA",
            alt="A",
            source_id="multi",
            vcf_source="panpop",
            alt_kind="sequence",
            allele_index=0,
            sample_genotypes={"GenomeA": "1/1", "GenomeB": "2/2"},
        )
        second_alt = SVRecord(
            genome="panpop",
            chrom="Chr1A",
            pos=100,
            end=100,
            svtype="INS",
            svlen=4,
            ref="A",
            alt="ATTTT",
            source_id="multi",
            vcf_source="panpop",
            alt_kind="sequence",
            allele_index=1,
            sample_genotypes={"GenomeA": "1/1", "GenomeB": "2/2"},
        )
        genome_a_pairwise = SVRecord(
            genome="GenomeA",
            chrom="Chr1A",
            pos=100,
            end=104,
            svtype="DEL",
            svlen=-5,
            ref="AAAAA",
            alt="A",
            vcf_source="pairwise",
            alt_kind="sequence",
        )
        genome_b_pairwise = SVRecord(
            genome="GenomeB",
            chrom="Chr1A",
            pos=100,
            end=100,
            svtype="INS",
            svlen=4,
            ref="A",
            alt="ATTTT",
            vcf_source="pairwise",
            alt_kind="sequence",
        )

        first_matches = find_matching_alt_records(
            first_alt,
            [genome_a_pairwise, genome_b_pairwise],
            criteria=MatchCriteria(max_distance=0, max_length_ratio_difference=0.0),
        )
        second_matches = find_matching_alt_records(
            second_alt,
            [genome_a_pairwise, genome_b_pairwise],
            criteria=MatchCriteria(max_distance=0, max_length_ratio_difference=0.0),
        )

        self.assertEqual(first_matches, [genome_a_pairwise])
        self.assertEqual(second_matches, [genome_b_pairwise])

    def test_extract_alt_sequence_can_use_query_assembly_flank_from_asm_info(self) -> None:
        panpop = SVRecord(
            genome="panpop",
            chrom="Chr1A",
            pos=100,
            end=100,
            svtype="INS",
            svlen=4,
            ref="A",
            alt="ATTTT",
            source_id="panpop-ins",
            vcf_source="panpop",
            alt_kind="sequence",
        )
        pairwise = SVRecord(
            genome="GenomeA",
            chrom="Chr1A",
            pos=100,
            end=100,
            svtype="INS",
            svlen=4,
            ref="A",
            alt="ATTTT",
            source_id="pairwise-ins",
            vcf_source="pairwise",
            alt_kind="sequence",
            info={
                "ASM_Chr": "1A",
                "ASM_Start": "4",
                "ASM_End": "7",
                "ASM_Strand": "+",
            },
        )
        query_sequences_by_genome = {"GenomeA": {"Chr1A": "AACCGGTTAA"}}

        extracted = extract_alt_sequence(
            panpop,
            pairwise,
            query_sequences_by_genome=query_sequences_by_genome,
            alt_flank=2,
        )

        self.assertEqual(extracted.sequence, "ACCGGTTA")
        self.assertEqual(extracted.md5, sequence_md5("ACCGGTTA"))

    def test_extract_alt_sequence_reverse_complements_negative_query_strand(self) -> None:
        panpop = SVRecord(
            genome="panpop",
            chrom="Chr1A",
            pos=100,
            end=100,
            svtype="INS",
            svlen=4,
            ref="A",
            alt="ATTTT",
            vcf_source="panpop",
            alt_kind="sequence",
        )
        pairwise = SVRecord(
            genome="GenomeA",
            chrom="Chr1A",
            pos=100,
            end=100,
            svtype="INS",
            svlen=4,
            ref="A",
            alt="ATTTT",
            vcf_source="pairwise",
            alt_kind="sequence",
            info={
                "ASM_Chr": "1A",
                "ASM_Start": "4",
                "ASM_End": "7",
                "ASM_Strand": "-",
            },
        )
        query_sequences_by_genome = {"GenomeA": {"1A": "AACCGGTTAA"}}

        extracted = extract_alt_sequence(
            panpop,
            pairwise,
            query_sequences_by_genome=query_sequences_by_genome,
            alt_flank=2,
        )

        self.assertEqual(extracted.sequence, "TAACCGGT")


    def test_min_panpop_sequence_length_filters_ref_and_alt_by_original_allele_length(self) -> None:
        panpop = SVRecord(
            genome="panpop",
            chrom="Chr1A",
            pos=2,
            end=2,
            svtype="INS",
            svlen=80,
            ref="A",
            alt="A" * 81,
            source_id="long-ins",
            vcf_source="panpop",
            alt_kind="sequence",
            sample_genotypes={"GenomeA": "1/1"},
        )
        pairwise = SVRecord(
            genome="GenomeA",
            chrom="Chr1A",
            pos=2,
            end=2,
            svtype="INS",
            svlen=80,
            ref="A",
            alt="C" * 81,
            source_id="pairwise-long-ins",
            vcf_source="pairwise",
            alt_kind="sequence",
        )

        extracted = list(
            iter_extracted_sequences(
                [panpop],
                [pairwise],
                {"Chr1A": "N" * 200},
                criteria=MatchCriteria(max_distance=0, max_length_ratio_difference=0.0),
                flank=50,
                min_panpop_sequence_length=50,
            )
        )

        self.assertEqual([record.role for record in extracted], ["alt"])
        self.assertEqual(extracted[0].sequence, "C" * 81)

    def test_min_panpop_sequence_length_filters_short_alt_but_keeps_long_ref(self) -> None:
        panpop = SVRecord(
            genome="panpop",
            chrom="Chr1A",
            pos=2,
            end=81,
            svtype="DEL",
            svlen=-80,
            ref="A" * 80,
            alt="A",
            source_id="long-del",
            vcf_source="panpop",
            alt_kind="sequence",
            sample_genotypes={"GenomeA": "1/1"},
        )
        pairwise = SVRecord(
            genome="GenomeA",
            chrom="Chr1A",
            pos=2,
            end=81,
            svtype="DEL",
            svlen=-80,
            ref="A" * 80,
            alt="A",
            source_id="pairwise-long-del",
            vcf_source="pairwise",
            alt_kind="sequence",
        )

        extracted = list(
            iter_extracted_sequences(
                [panpop],
                [pairwise],
                {"Chr1A": "N" * 200},
                criteria=MatchCriteria(max_distance=0, max_length_ratio_difference=0.0),
                flank=50,
                min_panpop_sequence_length=50,
            )
        )

        self.assertEqual([record.role for record in extracted], ["ref"])
        self.assertGreater(len(extracted[0].sequence), 80)
    def test_choose_representative_alt_uses_nearest_median_length(self) -> None:
        records = [
            ExtractedSequence(
                panpop_id="sv1",
                role="alt",
                genome="g1",
                chrom="Chr1",
                pos=1,
                end=1,
                svtype="INS",
                sequence="A" * 10,
                md5=sequence_md5("A" * 10),
            ),
            ExtractedSequence(
                panpop_id="sv1",
                role="alt",
                genome="g2",
                chrom="Chr1",
                pos=1,
                end=1,
                svtype="INS",
                sequence="A" * 20,
                md5=sequence_md5("A" * 20),
            ),
            ExtractedSequence(
                panpop_id="sv1",
                role="alt",
                genome="g3",
                chrom="Chr1",
                pos=1,
                end=1,
                svtype="INS",
                sequence="A" * 100,
                md5=sequence_md5("A" * 100),
            ),
        ]

        representative = choose_representative_alt(records)

        self.assertIsNotNone(representative)
        assert representative is not None
        self.assertEqual(representative.genome, "g2")
        self.assertEqual(representative.length, 20)


if __name__ == "__main__":
    unittest.main()

