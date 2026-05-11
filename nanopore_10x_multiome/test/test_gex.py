import pytest
from nanopore_10x_multiome.gex import (
    get_gex_anchors,
    TENX_GEX_ADAPTER,
    gex_re
)
from nanopore_10x_multiome.utils import (
    RC,
    get_barcode_parasail
)
from nanopore_10x_multiome.utils.test import (
    create_sequence,
    create_qual,
    BASE_SEQ
)

GEX_BARCODE = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCTATGCATGCATGCATGCATGCATGCATGCTTTTTTTTTTTTTTTTTTTTTT'
GBN = len(GEX_BARCODE)
BSN = len(BASE_SEQ)


def create_gex_sequence(
    barcode_seq,
    barcode_seq_loc,
    rc=False
):
    
    return create_sequence(
        barcode_seq,
        barcode_seq_loc,
        rc=rc
    )


def test_perfect_gex_parasail():
    seq = create_gex_sequence(
        GEX_BARCODE,
        0
    )

    psbc, psqual, psloc = get_barcode_parasail(
        seq,
        create_qual(len(seq), 33, 28),
        gex_re,
        TENX_GEX_ADAPTER,
        28
    )

    assert len(psbc) == 28
    assert psbc == "ATGCATGCATGCATGCATGCATGCATGC"
    assert psqual == "A" * 28
    assert psloc == 33


def test_extra_bases_gex_parasail():
    
    seq = create_gex_sequence(
        'ACACTCTTTCCCTTACACGGACGCTCCTTCCGATCTATGCATGCATGCATGCATGCATGCATGCTTTTTTTTTTTTTTTTTTTTTT',
        0
    )

    psbc, psqual, psloc = get_barcode_parasail(
        seq,
        create_qual(len(seq), 36, 28),
        gex_re,
        TENX_GEX_ADAPTER,
        28
    )

    assert len(psbc) == 28
    assert psbc == "ATGCATGCATGCATGCATGCATGCATGC"
    assert psqual == "A" * 28
    assert psloc == 36


def test_missing_bases_gex_parasail():
    seq = create_gex_sequence(
        'ACACTTTCCCTACACGACGCTCTTCCGATCTATGCATGCATGCATGCATGCATGCATGCTTTTTTTTTTTTTTTTTTTTTT',
        0
    )

    psbc, psqual, psloc = get_barcode_parasail(
        seq,
        create_qual(len(seq), 31, 28),
        gex_re,
        TENX_GEX_ADAPTER,
        28
    )

    assert len(psbc) == 28
    assert psbc == "ATGCATGCATGCATGCATGCATGCATGC"
    assert psqual == "A" * 28
    assert psloc == 31


def test_missing_bases_moved_gex_parasail():
    seq = create_gex_sequence(
        'ACACTCTTTCCCTACACGACGCTCTTCCGATCTATGCATGCATGCATGCATGCATGCATGCTTTTTTTTTTTTTTTTTTTTTT',
        550
    )

    psbc, psqual, psloc = get_barcode_parasail(
        seq,
        create_qual(len(seq), 550 + 33, 28),
        gex_re,
        TENX_GEX_ADAPTER,
        28
    )

    assert len(psbc) == 28
    assert psbc == "ATGCATGCATGCATGCATGCATGCATGC"
    assert psqual == "A" * 28
    assert psloc == 550 + 33


def test_perfect_gex_full():

    seq = create_gex_sequence(
        GEX_BARCODE,
        0
    )

    bc, umi, locs = get_gex_anchors(
        seq,
        create_qual(len(seq), 33, 16, 12)
    )

    genomic = seq[locs[0] + 22:locs[1]]

    assert locs == (GBN - 22, BSN + GBN)
    assert genomic == BASE_SEQ

    assert len(bc[0]) == 16
    assert bc[0] == "ATGCATGCATGCATGC"
    assert bc[1] == "A" * 16
    assert umi[0] == "ATGCATGCATGC"
    assert umi[1] == "B" * 12


def test_perfect_gex_full_rc():
    seq = create_gex_sequence(
        GEX_BARCODE,
        0,
        rc=True
    )

    bc, umi, locs = get_gex_anchors(
        seq,
        create_qual(len(seq), 33, 16, 12, rev=True),
    )

    genomic = seq[locs[0]:locs[1] - 22]

    assert genomic == RC(BASE_SEQ)
    assert locs == (0, BSN + 22)

    assert len(bc[0]) == 16
    assert bc[0] == "ATGCATGCATGCATGC"
    assert bc[1] == "A" * 16
    assert umi[0] == "ATGCATGCATGC"
    assert umi[1] == "B" * 12


def test_extra_bases_gex_full():

    seq = create_gex_sequence(
        'ACACTCTTTCCCTACACNGACGCTCNTTCCGATNCTATGCATGCATGCATGCATGCATGCATGCTTTTTTTTTTTTTTTTTTTTTT',
        0
    )

    bc, umi, locs = get_gex_anchors(
        seq,
        create_qual(len(seq), 36, 16, 12)
    )

    genomic = seq[locs[0] + 22:locs[1]]

    assert locs == (GBN - 22 + 3, BSN + GBN + 3)
    assert genomic == BASE_SEQ

    assert len(bc[0]) == 16
    assert bc[0] == "ATGCATGCATGCATGC"
    assert bc[1] == "A" * 16
    assert umi[0] == "ATGCATGCATGC"
    assert umi[1] == "B" * 12


def test_extra_bases_gex_full_rc():

    seq = create_gex_sequence(
        'ACACTCTTTCCCTACACNGACGCTCNTTCCGATNCTATGCATGCATGCATGCATGCATGCATGCTTTTTTTTTTTTTTTTTTTTTT',
        0,
        rc=True
    )

    bc, umi, locs = get_gex_anchors(
        seq,
        create_qual(len(seq), 36, 16, 12, rev=True)
    )

    genomic = seq[locs[0]:locs[1] - 22]

    assert locs == (0, BSN + 22)
    assert genomic == RC(BASE_SEQ)

    assert len(bc[0]) == 16
    assert bc[0] == "ATGCATGCATGCATGC"
    assert bc[1] == "A" * 16
    assert umi[0] == "ATGCATGCATGC"
    assert umi[1] == "B" * 12


def test_missing_bases_gex_full():

    seq = create_gex_sequence(
        'ACACTCTTTCCCTACACACGCTCTCCGATTATGCATGCATGCATGCATGCATGCATGCTTTTTTTTTTTTTTTTTTTTTT',
        0
    )

    bc, umi, locs = get_gex_anchors(
        seq,
        create_qual(len(seq), 30, 16, 12)
    )

    genomic = seq[locs[0] + 22:locs[1]]

    assert locs == (GBN - 22 - 3, BSN + GBN - 3)
    assert genomic == BASE_SEQ

    assert len(bc[0]) == 16
    assert bc[0] == "ATGCATGCATGCATGC"
    assert bc[1] == "A" * 16
    assert umi[0] == "ATGCATGCATGC"
    assert umi[1] == "B" * 12


def test_missing_bases_gex_full_rc_from_start():

    seq = create_gex_sequence(
        'ACACTCTTTCCCTACACACGCTCTCCGATTATGCATGCATGCATGCATGCATGCATGCTTTTTTTTTTTTTTTTTTTTTT',
        0,
        rc=True
    )

    bc, umi, locs = get_gex_anchors(
        seq,
        create_qual(len(seq), 30, 16, 12, rev=True)
    )

    genomic = seq[locs[0]:locs[1] - 22]

    assert locs == (0, BSN + 22)
    assert genomic == RC(BASE_SEQ)

    assert len(bc[0]) == 16
    assert bc[0] == "ATGCATGCATGCATGC"
    assert bc[1] == "A" * 16
    assert umi[0] == "ATGCATGCATGC"
    assert umi[1] == "B" * 12


def test_missing_bases_gex_full_at_start():

    seq = create_gex_sequence(
        'ACACTCTTTCCCTACACACGCTCTCCGATTATGCATGCATGCATGCATGCATGCATGCTTTTTTTTTTTTTTTTTTTTTT',
        0
    )

    bc, umi, locs = get_gex_anchors(
        seq,
        create_qual(len(seq), 30, 16, 12)
    )

    genomic = seq[locs[0] + 22:locs[1]]

    assert locs == (GBN - 22 - 3, BSN + GBN - 3)
    assert genomic == BASE_SEQ

    assert len(bc[0]) == 16
    assert bc[0] == "ATGCATGCATGCATGC"
    assert bc[1] == "A" * 16
    assert umi[0] == "ATGCATGCATGC"
    assert umi[1] == "B" * 12


def test_missing_bases_gex_full_rc():

    seq = create_gex_sequence(
        'ACACTCTTTCCCTACACACGCTCTCCGATTATGCATGCATGCATGCATGCATGCATGCTTTTTTTTTTTTTTTTTTTTTT',
        0,
        rc=True
    )

    bc, umi, locs = get_gex_anchors(
        seq,
        create_qual(len(seq), 30, 16, 12, rev=True)
    )

    genomic = seq[locs[0]:locs[1] - 22]

    assert locs == (0, BSN + 22)
    assert genomic == RC(BASE_SEQ)

    assert len(bc[0]) == 16
    assert bc[0] == "ATGCATGCATGCATGC"
    assert bc[1] == "A" * 16
    assert umi[0] == "ATGCATGCATGC"
    assert umi[1] == "B" * 12


def test_near_barcode_bases_gex_full():

    seq = create_gex_sequence(
        'ACACTCTTTCCCTACACGACGCTCTTCCGATATGCATGCATGCATGCATGCATGCATGCTTTATTTTTTTTTTTTTTTTTTT',
        0
    )

    bc, umi, locs = get_gex_anchors(
        seq,
        create_qual(len(seq), 31, 16, 12)
    )

    genomic = seq[locs[0] + 1 + 22:locs[1]]

    assert locs == (GBN - 22 - 2, BSN + GBN - 1)
    assert genomic == BASE_SEQ

    assert len(bc[0]) == 16
    assert bc[0] == "ATGCATGCATGCATGC"
    assert bc[1] == "A" * 16
    assert umi[0] == "ATGCATGCATGC"
    assert umi[1] == "B" * 12
