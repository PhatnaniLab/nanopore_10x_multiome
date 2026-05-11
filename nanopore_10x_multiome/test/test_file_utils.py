import pytest
import tempfile
import os

import pysam

from nanopore_10x_multiome.utils import (
    file_opener,
    get_file_writer,
    write_fastq_record,
    write_bam_record
)


class TestFileOpener:

    def test_open_fastq_for_writing(self, tmp_path):
        """file_opener opens .fastq files in text mode."""
        path = str(tmp_path / "test.fastq")
        with file_opener(path, mode='w') as fh:
            fh.write("@read1\nACGT\n+\nIIII\n")

        with open(path, 'r') as fh:
            content = fh.read()
        assert "@read1" in content

    def test_open_fastq_for_reading(self, tmp_path):
        """file_opener reads .fastq files."""
        path = str(tmp_path / "test.fastq")
        with open(path, 'w') as fh:
            fh.write("@read1\nACGT\n+\nIIII\n")

        with file_opener(path, mode='r') as fh:
            content = fh.read()
        assert "ACGT" in content

    def test_open_fastq_gz(self, tmp_path):
        """file_opener detects .fastq.gz and uses gzip."""
        import gzip
        path = str(tmp_path / "test.fastq.gz")
        with file_opener(path, mode='wt') as fh:
            fh.write("@read1\nACGT\n+\nIIII\n")

        with gzip.open(path, 'rt') as fh:
            content = fh.read()
        assert "ACGT" in content

    def test_open_bam_for_writing(self, tmp_path):
        """file_opener opens .bam files with pysam."""
        path = str(tmp_path / "test.bam")
        header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000, 'SN': 'chr1'}]}
        with file_opener(path, mode='w', header=header) as fh:
            assert isinstance(fh, pysam.AlignmentFile)

    def test_open_bam_adds_b_mode(self, tmp_path):
        """file_opener adds 'b' to mode for BAM files."""
        path = str(tmp_path / "test.bam")
        header = {'HD': {'VN': '1.0'}}
        # mode='w' should become 'wb' internally
        with file_opener(path, mode='w', header=header) as fh:
            assert isinstance(fh, pysam.AlignmentFile)

    def test_unknown_format_raises(self, tmp_path):
        """file_opener raises ValueError for unknown file extensions."""
        path = str(tmp_path / "test.xyz")
        with pytest.raises(ValueError, match="Unknown file format"):
            file_opener(path, mode='r')

    def test_explicit_file_format(self, tmp_path):
        """file_opener uses explicit file_format over extension."""
        path = str(tmp_path / "test.txt")
        with open(path, 'w') as fh:
            fh.write("@read1\nACGT\n+\nIIII\n")

        # Override with explicit format
        with file_opener(path, mode='r', file_format='fastq') as fh:
            content = fh.read()
        assert "ACGT" in content

    def test_explicit_gzip_flag(self, tmp_path):
        """file_opener respects gzip=True even without .gz extension."""
        import gzip
        path = str(tmp_path / "test.dat")
        with file_opener(path, mode='wt', file_format='fastq', gzip=True) as fh:
            fh.write("test data")

        with gzip.open(path, 'rt') as fh:
            assert fh.read() == "test data"


class TestGetFileWriter:

    def test_fastq_writer(self):
        """get_file_writer returns write_fastq_record for .fastq files."""
        writer = get_file_writer("output.fastq")
        assert writer is write_fastq_record

    def test_fastq_gz_writer(self):
        """get_file_writer returns write_fastq_record for .fastq.gz files."""
        writer = get_file_writer("output.fastq.gz")
        assert writer is write_fastq_record

    def test_bam_writer(self):
        """get_file_writer returns write_bam_record for .bam files."""
        writer = get_file_writer("output.bam")
        assert writer is write_bam_record

    def test_unknown_format_raises(self):
        """get_file_writer raises ValueError for unknown extensions."""
        with pytest.raises(ValueError, match="Unknown file format"):
            get_file_writer("output.xyz")

    def test_explicit_format_fastq(self):
        """get_file_writer uses explicit file_format."""
        writer = get_file_writer("anything.txt", file_format='fastq')
        assert writer is write_fastq_record

    def test_explicit_format_bam(self):
        """get_file_writer uses explicit file_format='bam'."""
        writer = get_file_writer("anything.txt", file_format='bam')
        assert writer is write_bam_record

    def test_explicit_unknown_format_raises(self):
        """get_file_writer raises ValueError for explicit unknown format."""
        with pytest.raises(ValueError, match="Unknown file format"):
            get_file_writer("output.fastq", file_format='csv')


class TestWriteFastqRecord:

    def test_basic_write(self, tmp_path):
        """write_fastq_record writes valid FASTQ format."""
        path = str(tmp_path / "test.fastq")
        with open(path, 'w') as fh:
            write_fastq_record(fh, "@read1", "ACGT", "IIII")

        with open(path, 'r') as fh:
            lines = fh.readlines()

        assert lines[0].strip() == "@read1"
        assert lines[1].strip() == "ACGT"
        assert lines[2].strip() == "+"
        assert lines[3].strip() == "IIII"

    def test_write_with_tags(self, tmp_path):
        """write_fastq_record appends tags to header."""
        path = str(tmp_path / "test.fastq")
        with open(path, 'w') as fh:
            write_fastq_record(fh, "@read1", "ACGT", "IIII", CB="AAAA", CR="AAAG")

        with open(path, 'r') as fh:
            header = fh.readline().strip()

        assert "CB=AAAA" in header
        assert "CR=AAAG" in header

    def test_none_tags_filtered(self, tmp_path):
        """write_fastq_record skips None-valued tags."""
        path = str(tmp_path / "test.fastq")
        with open(path, 'w') as fh:
            write_fastq_record(fh, "@read1", "ACGT", "IIII", CB="AAAA", CR=None)

        with open(path, 'r') as fh:
            header = fh.readline().strip()

        assert "CB=AAAA" in header
        assert "CR=" not in header

    def test_multiple_records(self, tmp_path):
        """Multiple records written sequentially."""
        path = str(tmp_path / "test.fastq")
        with open(path, 'w') as fh:
            write_fastq_record(fh, "@read1", "ACGT", "IIII")
            write_fastq_record(fh, "@read2", "TGCA", "FFFF")

        with open(path, 'r') as fh:
            lines = fh.readlines()

        assert len(lines) == 8  # 4 lines per record
        assert lines[4].strip() == "@read2"
