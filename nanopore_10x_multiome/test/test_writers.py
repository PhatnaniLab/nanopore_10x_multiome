import tempfile
import pysam
import os
from pathlib import Path

from nanopore_10x_multiome.utils import get_file_writer, fastqProcessor


TEST_FILE = os.path.join(Path(__file__).parent.absolute(), 'TEST_READS.fastq')


def test_fastq_writer():

    with tempfile.TemporaryDirectory() as td:

        _tempfile = os.path.join(td, 'temp.fastq')
        writer = get_file_writer(_tempfile)

        processor = fastqProcessor(
            verify_ids=False,
            phred_type='raw'
        )

        with open(TEST_FILE, mode='r') as fastqfile:
            with open(_tempfile, mode='w') as outfile:
                for x in processor.fastq_gen(fastqfile):
                    c, s, q = x[0]

                    writer(
                        outfile,
                        c,
                        s,
                        q
                    )

        with open(TEST_FILE, mode='r') as fastqfile:
            with open(_tempfile, mode='r') as outfile:
                for line1, line2 in zip(outfile, fastqfile):
                    assert line1 == line2


def test_bam_writer():

    with tempfile.TemporaryDirectory() as td:

        _tempfile = os.path.join(td, 'temp.bam')
        writer = get_file_writer(_tempfile)

        processor = fastqProcessor(
            verify_ids=False,
            phred_type='raw'
        )

        with open(TEST_FILE, mode='r') as fastqfile:
            with pysam.AlignmentFile(_tempfile, mode='wb', header={ 'HD': {'VN': '1.0'}}) as outfile:
                for x in processor.fastq_gen(fastqfile):
                    c, s, q = x[0]

                    writer(
                        outfile,
                        c,
                        s,
                        q
                    )

        processor = fastqProcessor(
            verify_ids=False
        )

        with open(TEST_FILE, mode='r') as fastqfile:
            with pysam.AlignmentFile(_tempfile, mode='rb', check_sq=False) as outfile:
                for x, y in zip(
                    processor.fastq_gen(fastqfile),
                    outfile
                ):
                    assert x[0][1] == y.query_sequence
                    assert x[0][2] == list(y.query_qualities)


def test_fastq_writer_tags():

    with tempfile.TemporaryDirectory() as td:

        _tempfile = os.path.join(td, 'temp.fastq')
        writer = get_file_writer(_tempfile)

        processor = fastqProcessor(
            verify_ids=False,
            phred_type='raw'
        )

        with open(TEST_FILE, mode='r') as fastqfile:
            with open(_tempfile, mode='w') as outfile:
                for x in processor.fastq_gen(fastqfile):
                    c, s, q = x[0]

                    writer(
                        outfile,
                        c,
                        s,
                        q,
                        **{'CB': 'ATGC', 'CR': 'ATTC', 'CY': 'II4I'}
                    )

        with open(TEST_FILE, mode='r') as fastqfile:
            with open(_tempfile, mode='r') as outfile:
                for line1, line2 in zip(processor.fastq_gen(outfile), processor.fastq_gen(fastqfile)):
                    assert line1[0][1] == line2[0][1]
                    assert line1[0][2] == line2[0][2]

                    assert ' CB=ATGC' in line1[0][0]
                    assert ' CR=ATTC' in line1[0][0]
                    assert ' CY=II4I' in line1[0][0]