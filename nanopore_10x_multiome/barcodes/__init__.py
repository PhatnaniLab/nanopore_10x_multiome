from ._load_barcodes import (
    load_atac_barcodes,
    load_gex_barcodes,
    load_translations,
    translate_barcode
)

from ._correct_barcodes import (
    correct_barcode,
    barcode_correction_table
)


def load_missing_multiome_barcode_info(pbar=False, test=False):
    BarcodeHolder.load(pbar=pbar, test=test)


class BarcodeHolder:

    gex_barcodes = None
    atac_barcodes = None
    gex_correction_table = None
    atac_correction_table = None
    atac_gex_translation_table = None

    @classmethod
    def load(cls, pbar=False, test=False):
        if test or cls.gex_barcodes is None:
            cls.gex_barcodes = load_gex_barcodes(test=test)

        if test or cls.atac_barcodes is None:
            cls.atac_barcodes = load_atac_barcodes(test=test)

        if cls.gex_correction_table is None:
            cls.gex_correction_table = barcode_correction_table(
                cls.gex_barcodes,
                pbar=pbar
            )

        if cls.atac_correction_table is None:
            cls.atac_correction_table = barcode_correction_table(
                cls.atac_barcodes,
                pbar=pbar
            )

        if cls.atac_gex_translation_table is None:
            cls.atac_gex_translation_table = load_translations(
                cls.atac_barcodes,
                cls.gex_barcodes
            )
