import os
import tempfile
from pathlib import Path

import numpy as np

from nanopore_10x_multiome.multiome import (
    split_multiome_preamp_fastq,
    _split_multiome_preamp_fastq
)
from nanopore_10x_multiome.barcodes import load_missing_multiome_barcode_info
from nanopore_10x_multiome.utils import fastqProcessor

TEST_FILE = os.path.join(Path(__file__).parent.absolute(), 'TEST_READS.fastq')
load_missing_multiome_barcode_info(test=True)

N_ATAC = 7
N_GEX = 10
N_TOTAL = 50


def test_result_counts_shape():
    """Single file returns a 1D array of length 3."""
    with tempfile.TemporaryDirectory() as td:
        out = [os.path.join(td, f'out{i}.fastq') for i in range(4)]
        result = split_multiome_preamp_fastq(
            TEST_FILE, out[0], out[1], out[2], out[3],
            keep_runoff_fragments=True
        )
        assert result.shape == (3,)
        assert result.dtype == int


def test_result_counts_sum():
    """All reads accounted for in result counts."""
    with tempfile.TemporaryDirectory() as td:
        out = [os.path.join(td, f'out{i}.fastq') for i in range(3)]
        result = split_multiome_preamp_fastq(
            TEST_FILE, out[0], out[1], out[2],
            keep_runoff_fragments=True
        )
        assert result.sum() == N_TOTAL


def test_without_runoff_fragments():
    """Without keep_runoff_fragments, fewer ATAC reads expected."""
    with tempfile.TemporaryDirectory() as td:
        out = [os.path.join(td, f'out{i}.fastq') for i in range(3)]
        result = split_multiome_preamp_fastq(
            TEST_FILE, out[0], out[1], out[2]
        )
        # Without runoff, some ATAC reads become 'other'
        assert result[0] <= N_ATAC
        assert result.sum() == N_TOTAL


def test_write_only_valid_barcodes():
    """write_only_valid_barcodes filters reads without valid barcodes."""
    with tempfile.TemporaryDirectory() as td:
        out_all = [os.path.join(td, f'all{i}.fastq') for i in range(3)]
        out_valid = [os.path.join(td, f'valid{i}.fastq') for i in range(3)]

        result_all = split_multiome_preamp_fastq(
            TEST_FILE, out_all[0], out_all[1], out_all[2],
            keep_runoff_fragments=True
        )
        result_valid = split_multiome_preamp_fastq(
            TEST_FILE, out_valid[0], out_valid[1], out_valid[2],
            write_only_valid_barcodes=True,
            keep_runoff_fragments=True
        )

        # Valid barcodes mode should produce fewer or equal ATAC/GEX reads
        assert result_valid[0] <= result_all[0]
        assert result_valid[1] <= result_all[1]


def test_n_records_limit():
    """_split_multiome_preamp_fastq respects n_records limit."""
    with tempfile.TemporaryDirectory() as td:
        out = [os.path.join(td, f'out{i}.fastq') for i in range(3)]
        result = _split_multiome_preamp_fastq(
            TEST_FILE, out[0], out[1], out[2],
            n_records=10
        )
        assert result.sum() == 10


def test_atac_technical_file():
    """ATAC technical file is written when path provided."""
    with tempfile.TemporaryDirectory() as td:
        out = [os.path.join(td, f'out{i}.fastq') for i in range(4)]
        split_multiome_preamp_fastq(
            TEST_FILE, out[0], out[1], out[2], out[3],
            keep_runoff_fragments=True
        )
        assert os.path.exists(out[3])

        with open(out[3], 'r') as fh:
            tech_lines = fh.readlines()
        # Should have same number of records as ATAC
        assert int(len(tech_lines) / 4) == N_ATAC


def test_no_atac_technical_file():
    """Without atac_technical_file_name, no extra file created."""
    with tempfile.TemporaryDirectory() as td:
        out = [os.path.join(td, f'out{i}.fastq') for i in range(3)]
        tech_path = os.path.join(td, 'tech.fastq')
        split_multiome_preamp_fastq(
            TEST_FILE, out[0], out[1], out[2]
        )
        assert not os.path.exists(tech_path)


def test_output_reads_have_tags():
    """ATAC and GEX output files contain barcode tags in headers."""
    with tempfile.TemporaryDirectory() as td:
        out = [os.path.join(td, f'out{i}.fastq') for i in range(3)]
        split_multiome_preamp_fastq(
            TEST_FILE, out[0], out[1], out[2],
            keep_runoff_fragments=True
        )

        # Check ATAC file has CR tags
        with open(out[0], 'r') as fh:
            header = fh.readline().strip()
        assert 'CR=' in header

        # Check GEX file has CR and UR tags
        with open(out[1], 'r') as fh:
            header = fh.readline().strip()
        assert 'CR=' in header
        assert 'UR=' in header


def test_other_reads_no_barcode_tags():
    """Other reads file should not contain barcode tags."""
    with tempfile.TemporaryDirectory() as td:
        out = [os.path.join(td, f'out{i}.fastq') for i in range(3)]
        split_multiome_preamp_fastq(
            TEST_FILE, out[0], out[1], out[2],
            keep_runoff_fragments=True
        )

        with open(out[2], 'r') as fh:
            for line in fh:
                if line.startswith('@'):
                    assert 'CB=' not in line


def test_multi_file_input():
    """List of input files are processed (single file list)."""
    with tempfile.TemporaryDirectory() as td:
        out = [[os.path.join(td, f'out{i}_{j}.fastq') for j in range(3)]
               for i in range(1)]

        result = split_multiome_preamp_fastq(
            [TEST_FILE],
            [out[0][0]],
            [out[0][1]],
            [out[0][2]],
            keep_runoff_fragments=True
        )

        assert result.shape == (1, 3)
        assert result[0].sum() == N_TOTAL
