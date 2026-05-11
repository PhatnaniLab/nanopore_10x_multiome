import pytest
import pysam

from nanopore_10x_multiome.utils._bam import bam_summarize_barcodes, write_bam_record


@pytest.fixture
def bam_with_barcodes(tmp_path):
    """Create a BAM file with reads having barcode tags."""
    bam_path = tmp_path / "test.bam"
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1000, 'SN': 'chr1'}]
    }

    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        write_bam_record(out, "read1", "ACGTACGT", "IIIIIIII", flag=0, CB="AAAA")
        write_bam_record(out, "read2", "ACGT", "IIII", flag=0, CB="BBBB")
        write_bam_record(out, "read3", "ACGTACGTACGT", "IIIIIIIIIIII", flag=0, CB="AAAA")
        write_bam_record(out, "read4", "ACGT", "IIII", flag=4, CB="CCCC")  # unmapped

    pysam.index(str(bam_path))
    return bam_path


def test_summarize_returns_dataframe(bam_with_barcodes):
    """bam_summarize_barcodes returns a DataFrame."""
    import pandas as pd
    result = bam_summarize_barcodes(bam_with_barcodes)
    assert isinstance(result, pd.DataFrame)


def test_summarize_columns(bam_with_barcodes):
    """Result has barcode, length, and is_mapped columns."""
    result = bam_summarize_barcodes(bam_with_barcodes)
    assert list(result.columns) == ['barcode', 'length', 'is_mapped']


def test_summarize_row_count(bam_with_barcodes):
    """One row per read in BAM file."""
    result = bam_summarize_barcodes(bam_with_barcodes)
    assert len(result) == 4


def test_summarize_barcodes(bam_with_barcodes):
    """Barcode values match what was written."""
    result = bam_summarize_barcodes(bam_with_barcodes)
    barcodes = list(result['barcode'])
    assert barcodes == ['AAAA', 'BBBB', 'AAAA', 'CCCC']


def test_summarize_lengths(bam_with_barcodes):
    """Lengths match sequence lengths."""
    result = bam_summarize_barcodes(bam_with_barcodes)
    lengths = list(result['length'])
    assert lengths == [8, 4, 12, 4]


def test_summarize_mapped_status(bam_with_barcodes):
    """Mapped status reflects alignment flags."""
    result = bam_summarize_barcodes(bam_with_barcodes)
    mapped = list(result['is_mapped'])
    # flag=4 means unmapped
    assert mapped[3] is False


def test_summarize_with_progress_bar(bam_with_barcodes):
    """Progress bar option works without changing results."""
    result_no_pbar = bam_summarize_barcodes(bam_with_barcodes, pbar=False)
    result_pbar = bam_summarize_barcodes(bam_with_barcodes, pbar=True)
    assert list(result_no_pbar['barcode']) == list(result_pbar['barcode'])


def test_summarize_custom_barcode_tag(tmp_path):
    """Custom barcode tag is used for extraction."""
    bam_path = tmp_path / "custom_tag.bam"
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1000, 'SN': 'chr1'}]
    }

    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        write_bam_record(out, "read1", "ACGT", "IIII", flag=0, BC="XXXX")

    pysam.index(str(bam_path))

    result = bam_summarize_barcodes(bam_path, barcode_tag='BC')
    assert result['barcode'].iloc[0] == "XXXX"


def test_summarize_empty_bam(tmp_path):
    """Empty BAM file produces empty DataFrame."""
    bam_path = tmp_path / "empty.bam"
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1000, 'SN': 'chr1'}]
    }

    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        pass

    pysam.index(str(bam_path))

    result = bam_summarize_barcodes(bam_path)
    assert len(result) == 0
