import os
from pathlib import Path
import tempfile

from nanopore_10x_multiome.multiome import split_multiome_preamp_fastq
from nanopore_10x_multiome.barcodes import load_missing_multiome_barcode_info
from nanopore_10x_multiome.utils import fastqProcessor

TEST_FILE = os.path.join(Path(__file__).parent.absolute(), 'TEST_READS.fastq')

N_ATAC = 7
N_GEX = 10

def test_multiome_stack():

    load_missing_multiome_barcode_info(test=True)

    with tempfile.TemporaryDirectory() as td:

        out_files = [
            os.path.join(td, f'out{i}.fastq')
            for i in range(4)
        ]

        split_multiome_preamp_fastq(
            TEST_FILE,
            out_files[0],
            out_files[1],
            out_files[2],
            out_files[3],
            keep_runoff_fragments=True
        )

        assert os.path.exists(out_files[0])
        assert os.path.exists(out_files[1])
        assert os.path.exists(out_files[2])
        assert os.path.exists(out_files[3])

        with open(out_files[1], mode='r') as test_file:
            assert N_GEX == int(len(list(test_file)) / 4)

        with open(out_files[0], mode='r') as test_file:
            assert N_ATAC == int(len(list(test_file)) / 4)

        with open(out_files[3], mode='r') as test_file:
            assert N_ATAC == int(len(list(test_file)) / 4)

        with open(out_files[2], mode='r') as test_file:
            assert (50 - N_GEX - N_ATAC) == int(len(list(test_file)) / 4)
