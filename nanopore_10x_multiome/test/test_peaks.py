import pytest
import pysam

from nanopore_10x_multiome.peaks.tn5_ends import (
    _trim_cigar_from_left,
    _trim_cigar_from_right,
    _make_single_end_read,
)

HEADER = pysam.AlignmentHeader.from_dict({
    'HD': {'VN': '1.6'},
    'SQ': [{'SN': 'chr1', 'LN': 1_000_000}],
})

# Op codes: M=0, I=1, D=2, N=3, S=4
M, I, D, N, S = 0, 1, 2, 3, 4


def _read(seq, cigar, ref_start, flag=0, name='r', qual=None, tags=None):
    r = pysam.AlignedSegment(HEADER)
    r.query_name = name
    r.query_sequence = seq
    r.query_qualities = pysam.qualitystring_to_array(
        qual if qual is not None else 'I' * len(seq)
    )
    r.flag = flag
    r.reference_id = 0
    r.reference_start = ref_start
    r.mapping_quality = 60
    r.cigar = cigar
    r.next_reference_id = 0
    r.next_reference_start = ref_start + 500
    r.template_length = 500
    if tags:
        r.tags = tags
    return r


# ── _trim_cigar_from_left ────────────────────────────────────────────────────

def test_trim_left_simple_match():
    cigar, ref = _trim_cigar_from_left([(M, 50)], 25)
    assert cigar == [(M, 25)]
    assert ref == 25


def test_trim_left_deletion_not_reached():
    # Window ends before the deletion
    cigar, ref = _trim_cigar_from_left([(M, 30), (D, 5), (M, 20)], 25)
    assert cigar == [(M, 25)]
    assert ref == 25


def test_trim_left_deletion_inside_window():
    # Window crosses the deletion: 30M + 5D + 5 of the trailing 20M
    cigar, ref = _trim_cigar_from_left([(M, 30), (D, 5), (M, 20)], 35)
    assert cigar == [(M, 30), (D, 5), (M, 5)]
    assert ref == 40


def test_trim_left_trailing_deletion_discarded():
    # The deletion falls exactly at the trailing edge and must be dropped
    cigar, ref = _trim_cigar_from_left([(M, 25), (D, 5)], 25)
    assert cigar == [(M, 25)]
    assert ref == 25


def test_trim_left_insertion_reduces_ref_consumed():
    # 10M5I35M  →  take 25 query bases = 10M + 5I + 10M; ref = 10 + 10 = 20
    cigar, ref = _trim_cigar_from_left([(M, 10), (I, 5), (M, 35)], 25)
    assert cigar == [(M, 10), (I, 5), (M, 10)]
    assert ref == 20


def test_trim_left_soft_clip_counted():
    # 5S45M  →  soft clip counts as query bases
    cigar, ref = _trim_cigar_from_left([(S, 5), (M, 45)], 25)
    assert cigar == [(S, 5), (M, 20)]
    assert ref == 20


# ── _trim_cigar_from_right ───────────────────────────────────────────────────

def test_trim_right_simple_match():
    cigar, ref = _trim_cigar_from_right([(M, 50)], 25)
    assert cigar == [(M, 25)]
    assert ref == 25


def test_trim_right_deletion_not_reached():
    # Window fits entirely in the rightmost block
    cigar, ref = _trim_cigar_from_right([(M, 30), (D, 5), (M, 20)], 10)
    assert cigar == [(M, 10)]
    assert ref == 10


def test_trim_right_deletion_inside_window():
    # BUG (current): produces "5M20M5D" instead of "5M5D20M"
    # Window of 25 spans: 5 from left block + 5D + all 20 of right block
    cigar, ref = _trim_cigar_from_right([(M, 30), (D, 5), (M, 20)], 25)
    assert cigar == [(M, 5), (D, 5), (M, 20)]
    assert ref == 30


def test_trim_right_leading_deletion_discarded():
    # Deletion is at the left (leading) edge of the original CIGAR;
    # it's past our window and must be dropped
    cigar, ref = _trim_cigar_from_right([(D, 5), (M, 30)], 20)
    assert cigar == [(M, 20)]
    assert ref == 20


def test_trim_right_insertion_reduces_ref_consumed():
    # 30M5I15M  →  last 10 query bases all come from the rightmost 15M block
    cigar, ref = _trim_cigar_from_right([(M, 30), (I, 5), (M, 15)], 10)
    assert cigar == [(M, 10)]
    assert ref == 10


def test_trim_right_insertion_spans_boundary():
    # 30M5I15M  →  last 25 = 5M + 5I + 15M; ref = 5 + 0 + 15 = 20
    # (insertion consumes query bases but not reference bases)
    cigar, ref = _trim_cigar_from_right([(M, 30), (I, 5), (M, 15)], 25)
    assert cigar == [(M, 5), (I, 5), (M, 15)]
    assert ref == 20


def test_trim_right_intron_inside_window():
    # 30M5D3N20M  →  last 25 = 5M + 5D + 3N + 20M; ref = 33
    cigar, ref = _trim_cigar_from_right([(M, 30), (D, 5), (N, 3), (M, 20)], 25)
    assert cigar == [(M, 5), (D, 5), (N, 3), (M, 20)]
    assert ref == 33


# ── _make_single_end_read ────────────────────────────────────────────────────

SEQ60 = 'A' * 25 + 'C' * 10 + 'T' * 25


def test_make_read_too_short_returns_none():
    r = _read('A' * 49, [(M, 49)], 100)
    assert _make_single_end_read(r, 25, HEADER) is None


def test_make_read_exact_minimum_length():
    # Exactly 2 * n_bp is the minimum that passes the length check
    r = _read('A' * 50, [(M, 50)], 100)
    assert _make_single_end_read(r, 25, HEADER) is not None


def test_make_read_returns_two_reads():
    r = _read(SEQ60, [(M, 60)], 100)
    outs = _make_single_end_read(r, 25, HEADER)
    assert outs is not None
    assert len(outs) == 2


def test_make_read_left_sequence():
    r = _read(SEQ60, [(M, 60)], 100)
    left, _ = _make_single_end_read(r, 25, HEADER)
    assert left.query_sequence == SEQ60[:25]


def test_make_read_right_sequence():
    r = _read(SEQ60, [(M, 60)], 100)
    _, right = _make_single_end_read(r, 25, HEADER)
    assert right.query_sequence == SEQ60[-25:]


def test_make_read_quality_scores_set():
    # BUG (current): query_qualities is never assigned, so it is None on output
    qual = 'F' * 25 + 'I' * 10 + '5' * 25  # 60 distinct quality chars
    r = _read(SEQ60, [(M, 60)], 100, qual=qual)
    left, right = _make_single_end_read(r, 25, HEADER)
    assert left.query_qualities is not None, "left read has no quality scores"
    assert right.query_qualities is not None, "right read has no quality scores"
    assert list(left.query_qualities) == list(pysam.qualitystring_to_array(qual[:25]))
    assert list(right.query_qualities) == list(pysam.qualitystring_to_array(qual[-25:]))


def test_make_read_left_reference_start():
    r = _read(SEQ60, [(M, 60)], 100)
    left, _ = _make_single_end_read(r, 25, HEADER)
    assert left.reference_start == 100


def test_make_read_right_reference_start_simple():
    r = _read(SEQ60, [(M, 60)], 100)
    _, right = _make_single_end_read(r, 25, HEADER)
    # reference_end = 100 + 60 = 160; last 25 bases → reference_start = 135
    assert right.reference_start == 135


def test_make_read_right_reference_start_with_deletion():
    # Read: 35M5D25M (query=60, ref=65), ref_start=100, ref_end=165
    seq = 'A' * 35 + 'T' * 25
    r = _read(seq, [(M, 35), (D, 5), (M, 25)], 100)
    left, right = _make_single_end_read(r, 25, HEADER)
    # Left: first 25 of 35M → ref_start=100
    assert left.reference_start == 100
    assert left.cigar == [(M, 25)]
    # Right: last 25 all within trailing 25M → ref_start = 165 - 25 = 140
    assert right.reference_start == 140
    assert right.cigar == [(M, 25)]


def test_make_read_paired_flags_cleared():
    flag = 0x1 | 0x2 | 0x8 | 0x20 | 0x40  # paired, proper, mate-unmapped, mate-rev, read1
    r = _read(SEQ60, [(M, 60)], 100, flag=flag)
    left, right = _make_single_end_read(r, 25, HEADER)
    for out in (left, right):
        assert not out.is_paired
        assert not out.is_proper_pair
        assert not out.is_read1
        assert not out.is_read2


def test_make_read_mate_fields_cleared():
    r = _read(SEQ60, [(M, 60)], 100)
    left, right = _make_single_end_read(r, 25, HEADER)
    for out in (left, right):
        assert out.next_reference_id == -1
        assert out.next_reference_start == -1
        assert out.template_length == 0


def test_make_read_mate_tags_stripped():
    r = _read(SEQ60, [(M, 60)], 100, tags=[('CB', 'ACGT'), ('MC', '60M'), ('MQ', 60)])
    left, right = _make_single_end_read(r, 25, HEADER)
    for out in (left, right):
        tag_keys = {t for t, _ in out.tags}
        assert 'CB' in tag_keys
        assert 'MC' not in tag_keys
        assert 'MQ' not in tag_keys


def test_make_read_query_name_preserved():
    r = _read(SEQ60, [(M, 60)], 100, name='cell1:frag42')
    left, right = _make_single_end_read(r, 25, HEADER)
    assert left.query_name == 'cell1:frag42'
    assert right.query_name == 'cell1:frag42'
