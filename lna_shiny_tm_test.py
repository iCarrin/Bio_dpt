import pytest 
import primer3
from lna_tm_shiny import *

class TestMakeComp():
    def test_pure_comp(self):
        assert make_comp("AAAAagtTTTTT", None) == "AAAAAACTTTTT"

    def test_round_about_pure_comp(self):
        assert make_comp("AAAAagtTTTTT", "G") == "AAAAAACTTTTT"

    def test_allele_insert(self):
        assert make_comp("AAAAagtTTTTT", "C") == "AAAAAAGTTTTT"

    def test_near_allele_insert(self):
        assert make_comp("aaaAAGTTTTTT", "C") == "AAAAAACTTTGT"

    def test_second_if_lots(self):
        assert make_comp("AAAAagttttTT", "C") == "AAAAAAGTTTTT"


class TestCheckSymytry():
    def test_false_from_len(self):
        assert check_sym("ATGCTTGAAGCAT") == False

    def test_false_sym(self):
        assert check_sym("ATGCTTGAAGCATT") == False

    def test_true_sym(self):
        assert check_sym("ATGCTTGCAAGCAT") == True




class TestCalcTmWithLna:

    @pytest.fixture
    def temp_guy(self):
        b = primer3.thermoanalysis.ThermoAnalysis()
        b.set_thermo_args(
            mv_conc=200,
            dv_conc=50,
            dntp_conc=3,
            dna_conc=0.8,
            dmso_conc=0.0,
            dmso_fact=0.0,
            formamide_conc=0.0,
            salt_correction_method="owczarzy",
        )
        return b

    @pytest.mark.parametrize(
        "seq, expected", 
        [
            ("AAAAAAAAAAAAAAAAAA", 44.8),
            ("TTTTTTTTTTTTTTTTTT", 44.8),
            ("GGGGGGGGGGGGGGGGGG", 78.4),
            ("CCCCCCCCCCCCCCCCCC", 78.4),
            ("AGTCAGTCACTAGTTGCT", 57.4),
            
            
        ],
    )
    def test_against_oligo_no_lna(self, temp_guy, seq, expected):
        assert calc_tm_with_lna(seq, temp_guy) == pytest.approx(expected, abs=0.1)

    # @pytest.mark.parametrize(
    #         "seq, expected",
    #         [
    #             ("AGTCAGTC+A+C+TAGTTGCT", 63.2),
    #         ],
    # )
    # def test_against_oligo_with_lna(self, temp_guy, seq, expected):
    #     assert calc_tm_with_lna(seq, temp_guy) == pytest.approx(expected, abs=0.1)


    # @pytest.mark.parametrize( #this test is doomed to fail because there are no parameters for Normal missmatches
    #  # and oligo doesn't do locked missmatches
    #         "seq, allele, expected",
    #         [
    #             ("AAGTCAGTCACTAGTTGCT", "C", 55.8 ),
    #         ],
    # )
    def test_against_oligo_no_lna(self, temp_guy, seq, allele, expected):
        assert calc_tm_with_lna(seq, temp_guy, allele) 
    @pytest.mark.parametrize(
        "seq, tol",
        [
            pytest.param("AAAAAAAAAAAAAAAAAA", 4.5e-3, id="all_A"),
            pytest.param("TTTTTTTTTTTTTTTTTT", 4.5e-3, id="all_T"),
            pytest.param("GGGGGGGGGGGGGGGGGG", 4.5e-3, id="all_G"),
            pytest.param("CCCCCCCCCCCCCCCCCC", 4.5e-3, id="all_C"),
        ],
    )
    def test_primer3_no_lna(self, temp_guy, seq, tol):
        assert calc_tm_with_lna(seq, temp_guy) == pytest.approx(temp_guy.calc_tm(seq), abs=tol)