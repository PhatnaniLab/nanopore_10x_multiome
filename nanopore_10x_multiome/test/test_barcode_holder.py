import pytest
import numpy as np

from nanopore_10x_multiome.barcodes import (
    BarcodeHolder,
    load_missing_multiome_barcode_info
)


@pytest.fixture(autouse=True)
def reset_barcode_holder():
    """Reset BarcodeHolder state before each test."""
    BarcodeHolder.gex_barcodes = None
    BarcodeHolder.atac_barcodes = None
    BarcodeHolder.gex_correction_table = None
    BarcodeHolder.atac_correction_table = None
    BarcodeHolder.atac_gex_translation_table = None
    yield
    # Reset again after test
    BarcodeHolder.gex_barcodes = None
    BarcodeHolder.atac_barcodes = None
    BarcodeHolder.gex_correction_table = None
    BarcodeHolder.atac_correction_table = None
    BarcodeHolder.atac_gex_translation_table = None


def test_load_populates_all_fields():
    """Loading BarcodeHolder fills all class attributes."""
    BarcodeHolder.load(test=True)

    assert BarcodeHolder.gex_barcodes is not None
    assert BarcodeHolder.atac_barcodes is not None
    assert BarcodeHolder.gex_correction_table is not None
    assert BarcodeHolder.atac_correction_table is not None
    assert BarcodeHolder.atac_gex_translation_table is not None


def test_load_test_mode_barcode_count():
    """Test mode loads 100 barcodes for each type."""
    BarcodeHolder.load(test=True)

    assert len(BarcodeHolder.gex_barcodes) == 100
    assert len(BarcodeHolder.atac_barcodes) == 100


def test_singleton_caching():
    """Second load call does not reload (cached)."""
    BarcodeHolder.load(test=True)
    gex_barcodes_first = BarcodeHolder.gex_barcodes

    BarcodeHolder.load()
    # Same object - not reloaded
    assert BarcodeHolder.gex_barcodes is gex_barcodes_first


def test_correction_tables_are_dicts():
    """Correction tables are dictionaries."""
    BarcodeHolder.load(test=True)

    assert isinstance(BarcodeHolder.gex_correction_table, dict)
    assert isinstance(BarcodeHolder.atac_correction_table, dict)


def test_correction_tables_contain_originals():
    """Correction tables contain all original barcodes mapping to themselves."""
    BarcodeHolder.load(test=True)

    for bc in BarcodeHolder.gex_barcodes:
        assert BarcodeHolder.gex_correction_table[bc] == bc

    for bc in BarcodeHolder.atac_barcodes:
        assert BarcodeHolder.atac_correction_table[bc] == bc


def test_translation_table_size():
    """Translation table maps all ATAC barcodes to GEX barcodes."""
    BarcodeHolder.load(test=True)

    assert len(BarcodeHolder.atac_gex_translation_table) == 100


def test_translation_table_maps_atac_to_gex():
    """Translation table values are GEX barcodes."""
    BarcodeHolder.load(test=True)

    gex_set = set(BarcodeHolder.gex_barcodes)
    for atac_bc, gex_bc in BarcodeHolder.atac_gex_translation_table.items():
        assert atac_bc in BarcodeHolder.atac_barcodes
        assert gex_bc in gex_set


def test_load_missing_wrapper():
    """load_missing_multiome_barcode_info wrapper works."""
    load_missing_multiome_barcode_info(test=True)

    assert BarcodeHolder.gex_barcodes is not None
    assert BarcodeHolder.atac_barcodes is not None


def test_partial_load_fills_missing():
    """If only some fields are loaded, load fills the rest."""
    from nanopore_10x_multiome.barcodes._load_barcodes import load_gex_barcodes
    BarcodeHolder.gex_barcodes = load_gex_barcodes(test=True)

    BarcodeHolder.load(test=True)

    # All fields should be populated
    assert BarcodeHolder.atac_barcodes is not None
    assert BarcodeHolder.gex_correction_table is not None
    assert BarcodeHolder.atac_correction_table is not None
    assert BarcodeHolder.atac_gex_translation_table is not None
