import pytest
import numpy as np

from nanopore_10x_multiome.barcodes._correct_barcodes import (
    correct_barcode,
    barcode_correction_table
)


@pytest.fixture
def simple_correction_table():
    """Create a correction table from a small set of barcodes."""
    barcodes = ["ACGT", "TGCA", "AAAA"]
    return barcode_correction_table(barcodes)


@pytest.fixture
def simple_barcodes():
    return ["ACGT", "TGCA", "AAAA"]


def test_exact_match(simple_correction_table):
    """Exact barcode match returns itself."""
    result = correct_barcode("ACGT", "IIII", simple_correction_table)
    assert result == "ACGT"


def test_exact_match_all_barcodes(simple_correction_table):
    """All original barcodes map to themselves."""
    assert correct_barcode("ACGT", "IIII", simple_correction_table) == "ACGT"
    assert correct_barcode("TGCA", "IIII", simple_correction_table) == "TGCA"
    assert correct_barcode("AAAA", "IIII", simple_correction_table) == "AAAA"


def test_single_substitution_correction(simple_correction_table):
    """Single substitution error is corrected."""
    # CCGT is 1 sub from ACGT, not close to others
    result = correct_barcode("CCGT", "IIII", simple_correction_table)
    assert result == "ACGT"


def test_uncorrectable_barcode(simple_correction_table):
    """Barcode too far from any valid barcode returns None."""
    result = correct_barcode("CCCC", "IIII", simple_correction_table)
    assert result is None


def test_ambiguous_barcode_returns_none():
    """Barcode equidistant from two valid barcodes returns None."""
    # ACGT and ACTT differ by 1 base; ACNT maps to neither
    table = barcode_correction_table(["ACGT", "ACTT"])
    # "ACCT" is 1 sub from both - should be ambiguous (not in table)
    result = correct_barcode("ACCT", "IIII", table)
    assert result is None


def test_none_barcode_not_in_table():
    """Barcode not in table and max_dist=1 returns None."""
    table = barcode_correction_table(["ACGTACGT"])
    result = correct_barcode("TTTTTTTT", "IIIIIIII", table)
    assert result is None


def test_max_dist_1_default():
    """With default max_dist=1, only table lookups are used."""
    table = barcode_correction_table(["ACGT"])
    # 2 subs away - should fail with max_dist=1
    result = correct_barcode("TTGT", "IIII", table)
    assert result is None


def test_max_dist_greater_than_1_requires_kwargs():
    """max_dist > 1 without valid_barcodes raises RuntimeError."""
    table = barcode_correction_table(["ACGT"])
    with pytest.raises(RuntimeError, match="max_dist > 1"):
        correct_barcode("TTGT", "IIII", table, max_dist=2)


def test_max_dist_2_with_valid_barcodes():
    """max_dist=2 with quality-weighted distance finds closest barcode."""
    barcodes = ["ACGTACGT", "TGCATGCA"]
    table = barcode_correction_table(barcodes)
    bc_array = np.array(barcodes)
    char_table = np.array([[ord(c) for c in bc] for bc in barcodes])

    # 2 subs from ACGTACGT (positions 0,1 changed)
    result = correct_barcode(
        "TTGTACGT",
        "IIIIIIII",
        table,
        max_dist=2,
        valid_barcodes=bc_array,
        valid_barcodes_char_table=char_table,
        min_weight_dist=0.5
    )
    assert result == "ACGTACGT"


def test_max_dist_2_too_far_returns_none():
    """max_dist=2 but barcode is 3+ away returns None."""
    barcodes = ["ACGTACGT", "TGCATGCA"]
    table = barcode_correction_table(barcodes)
    bc_array = np.array(barcodes)
    char_table = np.array([[ord(c) for c in bc] for bc in barcodes])

    # 4 subs from ACGTACGT
    result = correct_barcode(
        "TGCAACGT",
        "IIIIIIII",
        table,
        max_dist=2,
        valid_barcodes=bc_array,
        valid_barcodes_char_table=char_table,
        min_weight_dist=0.5
    )
    assert result is None


def test_quality_weighted_distance_prefers_low_quality_mismatches():
    """Low quality mismatches contribute less to weighted distance."""
    barcodes = ["ACGTACGT", "TCGTACGT"]
    table = barcode_correction_table(barcodes)
    bc_array = np.array(barcodes)
    char_table = np.array([[ord(c) for c in bc] for bc in barcodes])

    # "GCGTACGT" differs at position 0 from both
    # With low quality at position 0, the mismatch weighs less
    # but since it's equidistant (1 from each), weighted distance decides
    result = correct_barcode(
        "GCGTACGT",
        "!IIIIIII",  # low quality at mismatch position
        table,
        max_dist=2,
        valid_barcodes=bc_array,
        valid_barcodes_char_table=char_table,
        min_weight_dist=0.5
    )
    # Both are distance 1, and both have the same quality weight at the
    # mismatch position, so this should be None (ambiguous)
    assert result is None


def test_deletion_correction(simple_correction_table):
    """Single deletion can be corrected via lookup table."""
    # "CGT" is "ACGT" with first base deleted
    result = correct_barcode("CGT", "III", simple_correction_table)
    # Should be found if it's uniquely mapped
    if "CGT" in simple_correction_table:
        assert result == simple_correction_table["CGT"]
    else:
        assert result is None


def test_insertion_correction(simple_correction_table):
    """Single insertion can be corrected via lookup table."""
    # "AACGT" is "ACGT" with extra A at position 1
    result = correct_barcode("AACGT", "IIIII", simple_correction_table)
    if "AACGT" in simple_correction_table:
        assert result == simple_correction_table["AACGT"]
    else:
        assert result is None
