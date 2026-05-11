import pytest

from nanopore_10x_multiome.atac import process_atac_tags
from nanopore_10x_multiome.gex import process_gex_tags
from nanopore_10x_multiome.barcodes._correct_barcodes import barcode_correction_table
from nanopore_10x_multiome.barcodes._load_barcodes import load_translations


# ATAC tag processing tests

class TestProcessAtacTags:

    @pytest.fixture
    def atac_setup(self):
        """Set up ATAC barcodes and correction/translation tables."""
        atac_barcodes = ["ACGTACGTACGTACGT", "TGCATGCATGCATGCA"]
        gex_barcodes = ["AAAACCCCGGGGTTTT", "TTTTGGGGCCCCAAAA"]
        correction_table = barcode_correction_table(atac_barcodes)
        translation_table = dict(zip(atac_barcodes, gex_barcodes))
        return correction_table, translation_table

    def test_valid_barcode(self, atac_setup):
        """Valid ATAC barcode produces correct tags and is_valid=True."""
        correction_table, translation_table = atac_setup
        tags, is_valid = process_atac_tags(
            "ACGTACGTACGTACGT",
            "IIIIIIIIIIIIIIII",
            correction_table,
            translation_table
        )

        assert is_valid is True
        assert tags['CB'] == "AAAACCCCGGGGTTTT"  # translated to GEX
        assert tags['CR'] == "ACGTACGTACGTACGT"  # raw barcode preserved
        assert tags['CY'] == "IIIIIIIIIIIIIIII"  # quality preserved

    def test_correctable_barcode(self, atac_setup):
        """Single-error barcode is corrected and translated."""
        correction_table, translation_table = atac_setup
        # Single sub at position 0: TCGTACGTACGTACGT
        tags, is_valid = process_atac_tags(
            "TCGTACGTACGTACGT",
            "IIIIIIIIIIIIIIII",
            correction_table,
            translation_table
        )

        assert is_valid is True
        assert tags['CB'] == "AAAACCCCGGGGTTTT"
        assert tags['CR'] == "TCGTACGTACGTACGT"

    def test_uncorrectable_barcode(self, atac_setup):
        """Barcode too far from any valid returns CB=None, is_valid=False."""
        correction_table, translation_table = atac_setup
        tags, is_valid = process_atac_tags(
            "GGGGGGGGGGGGGGGG",
            "IIIIIIIIIIIIIIII",
            correction_table,
            translation_table
        )

        assert is_valid is False
        assert tags['CB'] is None
        assert tags['CR'] == "GGGGGGGGGGGGGGGG"
        assert tags['CY'] == "IIIIIIIIIIIIIIII"

    def test_tags_keys(self, atac_setup):
        """ATAC tags always include CB, CR, CY."""
        correction_table, translation_table = atac_setup
        tags, _ = process_atac_tags(
            "ACGTACGTACGTACGT",
            "IIIIIIIIIIIIIIII",
            correction_table,
            translation_table
        )

        assert set(tags.keys()) == {'CB', 'CR', 'CY'}

    def test_barcode_not_in_translation_table(self):
        """Corrected barcode not in translation table passes through."""
        barcode = "ACGTACGTACGTACGT"
        correction_table = {barcode: barcode}
        translation_table = {}  # empty - no translations

        tags, is_valid = process_atac_tags(
            barcode,
            "IIIIIIIIIIIIIIII",
            correction_table,
            translation_table
        )

        assert is_valid is True
        # translate_barcode returns the barcode itself if not in table
        assert tags['CB'] == barcode


# GEX tag processing tests

class TestProcessGexTags:

    @pytest.fixture
    def gex_setup(self):
        """Set up GEX barcodes and correction table."""
        gex_barcodes = ["ACGTACGTACGTACGT", "TGCATGCATGCATGCA"]
        correction_table = barcode_correction_table(gex_barcodes)
        return correction_table

    def test_valid_barcode_and_umi(self, gex_setup):
        """Valid GEX barcode returns correct tags with UMI."""
        correction_table = gex_setup
        tags, is_valid = process_gex_tags(
            "ACGTACGTACGTACGT",
            "IIIIIIIIIIIIIIII",
            "ATGCATGCATGC",
            "FFFFFFFFFFFF",
            correction_table
        )

        assert is_valid is True
        assert tags['CB'] == "ACGTACGTACGTACGT"
        assert tags['CR'] == "ACGTACGTACGTACGT"
        assert tags['CY'] == "IIIIIIIIIIIIIIII"
        assert tags['UB'] == "ATGCATGCATGC"
        assert tags['UR'] == "ATGCATGCATGC"
        assert tags['UY'] == "FFFFFFFFFFFF"

    def test_correctable_barcode(self, gex_setup):
        """Single-error GEX barcode is corrected."""
        correction_table = gex_setup
        tags, is_valid = process_gex_tags(
            "TCGTACGTACGTACGT",
            "IIIIIIIIIIIIIIII",
            "ATGCATGCATGC",
            "FFFFFFFFFFFF",
            correction_table
        )

        assert is_valid is True
        assert tags['CB'] == "ACGTACGTACGTACGT"
        assert tags['CR'] == "TCGTACGTACGTACGT"

    def test_uncorrectable_barcode(self, gex_setup):
        """Uncorrectable GEX barcode returns is_valid=False."""
        correction_table = gex_setup
        tags, is_valid = process_gex_tags(
            "GGGGGGGGGGGGGGGG",
            "IIIIIIIIIIIIIIII",
            "ATGCATGCATGC",
            "FFFFFFFFFFFF",
            correction_table
        )

        assert is_valid is False
        assert tags['CB'] is None
        assert tags['CR'] == "GGGGGGGGGGGGGGGG"

    def test_tags_keys(self, gex_setup):
        """GEX tags always include CB, CR, CY, UB, UR, UY."""
        correction_table = gex_setup
        tags, _ = process_gex_tags(
            "ACGTACGTACGTACGT",
            "IIIIIIIIIIIIIIII",
            "ATGCATGCATGC",
            "FFFFFFFFFFFF",
            correction_table
        )

        assert set(tags.keys()) == {'CB', 'CR', 'CY', 'UB', 'UR', 'UY'}

    def test_umi_passthrough(self, gex_setup):
        """UMI values are passed through without correction."""
        correction_table = gex_setup
        tags, _ = process_gex_tags(
            "ACGTACGTACGTACGT",
            "IIIIIIIIIIIIIIII",
            "NNNNNNNNNNN",
            "!!!!!!!!!!!",
            correction_table
        )

        # UB and UR should both be the raw UMI
        assert tags['UB'] == "NNNNNNNNNNN"
        assert tags['UR'] == "NNNNNNNNNNN"
        assert tags['UY'] == "!!!!!!!!!!!"
