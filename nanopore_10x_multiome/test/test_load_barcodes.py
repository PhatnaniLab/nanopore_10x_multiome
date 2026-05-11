import pytest
import numpy as np

from nanopore_10x_multiome.barcodes._load_barcodes import (
    load_atac_barcodes,
    load_gex_barcodes,
    load_translations,
    translate_barcode
)


def test_load_atac_barcodes_test_mode():
    """Loading ATAC barcodes in test mode returns exactly 100."""
    barcodes = load_atac_barcodes(test=True)
    assert len(barcodes) == 100
    assert isinstance(barcodes, np.ndarray)


def test_load_gex_barcodes_test_mode():
    """Loading GEX barcodes in test mode returns exactly 100."""
    barcodes = load_gex_barcodes(test=True)
    assert len(barcodes) == 100
    assert isinstance(barcodes, np.ndarray)


def test_atac_barcodes_are_strings():
    """ATAC barcodes should be DNA strings."""
    barcodes = load_atac_barcodes(test=True)
    for bc in barcodes:
        assert isinstance(bc, str)
        assert all(c in 'ACGTN' for c in bc.upper())


def test_gex_barcodes_are_strings():
    """GEX barcodes should be DNA strings."""
    barcodes = load_gex_barcodes(test=True)
    for bc in barcodes:
        assert isinstance(bc, str)
        assert all(c in 'ACGTN' for c in bc.upper())


def test_atac_and_gex_barcodes_differ():
    """ATAC and GEX barcodes are distinct sets."""
    atac = set(load_atac_barcodes(test=True))
    gex = set(load_gex_barcodes(test=True))
    assert atac != gex


def test_atac_barcodes_uniform_length():
    """All ATAC barcodes have the same length (16bp)."""
    barcodes = load_atac_barcodes(test=True)
    lengths = set(len(bc) for bc in barcodes)
    assert len(lengths) == 1
    assert lengths.pop() == 16


def test_gex_barcodes_uniform_length():
    """All GEX barcodes have the same length (16bp)."""
    barcodes = load_gex_barcodes(test=True)
    lengths = set(len(bc) for bc in barcodes)
    assert len(lengths) == 1
    assert lengths.pop() == 16


def test_load_translations_with_barcodes():
    """Translation table maps ATAC barcodes to GEX barcodes positionally."""
    atac = load_atac_barcodes(test=True)
    gex = load_gex_barcodes(test=True)
    table = load_translations(atac, gex)

    assert len(table) == 100
    # First ATAC barcode maps to first GEX barcode
    assert table[atac[0]] == gex[0]
    # Last ATAC barcode maps to last GEX barcode
    assert table[atac[-1]] == gex[-1]


def test_load_translations_default_loading():
    """Translation table loads barcodes from files when not provided."""
    table = load_translations()
    assert len(table) > 0
    # All values should be strings
    for k, v in list(table.items())[:10]:
        assert isinstance(k, str)
        assert isinstance(v, str)


def test_translate_barcode_valid():
    """Valid ATAC barcode translates to GEX barcode."""
    atac = load_atac_barcodes(test=True)
    gex = load_gex_barcodes(test=True)
    table = load_translations(atac, gex)

    result = translate_barcode(atac[0], table)
    assert result == gex[0]


def test_translate_barcode_none():
    """None barcode returns None."""
    table = {"ACGT": "TGCA"}
    result = translate_barcode(None, table)
    assert result is None


def test_translate_barcode_not_in_table():
    """Unknown barcode returns itself (passthrough)."""
    table = {"ACGT": "TGCA"}
    result = translate_barcode("GGGG", table)
    assert result == "GGGG"


def test_translate_barcode_all_in_table():
    """Every ATAC barcode in test set can be translated."""
    atac = load_atac_barcodes(test=True)
    gex = load_gex_barcodes(test=True)
    table = load_translations(atac, gex)

    for a, g in zip(atac, gex):
        assert translate_barcode(a, table) == g
