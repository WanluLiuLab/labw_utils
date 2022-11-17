"""
test_fasta.py -- Unit test of corresponding module.
"""
import itertools
import os
import tempfile

import pytest

from labw_utils.bioutils.datastructure.fai_view import create_fai_from_fasta
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory, FastaViewType, SeekTooFarError, \
    ChromosomeNotFoundError, \
    FromGreaterThanToError, DuplicatedChromosomeNameError
from labw_utils.bioutils.parser.fai import FastaIndexNotWritableError
from labw_utils.commonutils.io.safe_io import get_writer, get_reader
from labw_utils.commonutils.stdlib_helper import logger_helper
from test_labw_utils.bioutils import TEST_DATA_DIR

lh = logger_helper.get_logger(__name__)

test_fasta_path = os.path.join(TEST_DATA_DIR, "test.fasta")
test_oneline_fasta_path = os.path.join(TEST_DATA_DIR, "test_oneline.fasta")

with get_reader(test_fasta_path) as reader:
    test_fasta_seq = reader.read()

with get_reader(test_oneline_fasta_path) as reader:
    test_oneline_fasta_seq = reader.read()

VALID_NEWLINE = ("\n", "\r\n")


def _public_asserts(fa: FastaViewType) -> None:
    assert len(fa) == 5
    assert fa.sequence('chr3', 2, 15) == 'TANNTGNATNATG'  # At line end
    assert fa.sequence('chr3', 2, 16) == 'TANNTGNATNATGN'  # Cross the line
    assert fa.sequence('chr2') == 'NNNNNNNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTTCTNNNNNNNNN'
    assert fa.sequence('chr2') == fa.sequence('chr2', 0)
    assert fa.sequence('chr2') == fa.sequence('chr2', 0, -1)
    assert fa.sequence('chr4') == 'AAAAAAAAAACCCCCC'
    assert fa.sequence('chr6') == 'CTA'
    assert fa.query(('chr3', 2, 15)) == fa.sequence('chr3', 2, 15)
    assert fa.query(('chr3', 2)) == fa.sequence('chr3', 2)
    assert fa.query(('chr3',), ) == fa.sequence('chr3')

    with pytest.raises(SeekTooFarError):
        fa.sequence('chr2', 5, 1222)
    with pytest.raises(SeekTooFarError):
        fa.sequence('chr2', -5, 29)
    with pytest.raises(SeekTooFarError):
        fa.sequence('chr2', 500, 29)
    with pytest.raises(FromGreaterThanToError):
        fa.sequence('chr2', 12, 5)


def _fasta_with_full_header_assets(fa: FastaViewType) -> None:
    assert fa.sequence('chr1 some att', 0, 1) == 'N'
    assert fa.sequence('chr1 some att', 26, 29) == 'CCA'  # Cross the line
    assert fa.sequence('chr1 some att', 28, 29) == 'A'  # Next line
    assert fa.sequence('chr1 some att', 5, 29) == 'NNNNNNNNNNATCGTTACGTACCA'
    with pytest.raises(ChromosomeNotFoundError):
        fa.sequence('chr1', 5, 29)
    assert fa.sequence('chr1 some att', 5, 63) == 'NNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTT'
    _public_asserts(fa)


def _fasta_without_full_header_assets(fa: FastaViewType) -> None:
    assert fa.sequence('chr1', 0, 1) == 'N'
    assert fa.sequence('chr1', 26, 29) == 'CCA'  # Cross the line
    assert fa.sequence('chr1', 28, 29) == 'A'  # Next line
    assert fa.sequence('chr1', 5, 29) == 'NNNNNNNNNNATCGTTACGTACCA'
    with pytest.raises(ChromosomeNotFoundError):
        fa.sequence('chr1 some att', 5, 29)
    assert fa.sequence('chr1', 5, 63) == 'NNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTT'
    _public_asserts(fa)


@pytest.mark.parametrize(
    argnames="kwargs",
    argvalues=(
            {
                i: j for i, j in zip(
                ["read_into_memory", "full_header", "newline", "fasta_seq"],
                test_case
            )
            } for test_case in itertools.product(
        [True, False],
        [True, False],
        VALID_NEWLINE,
        [test_oneline_fasta_seq, test_fasta_seq]
    )
    )
)
def test_newline_with_or_without_full_header(kwargs) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_filename = os.path.join(tmpdir, "tmp.fa")
        with get_writer(fasta_filename, newline=kwargs["newline"]) as fh:
            fh.write(kwargs["fasta_seq"])
        with FastaViewFactory(
                filename=fasta_filename,
                full_header=kwargs['full_header'],
                read_into_memory=kwargs["read_into_memory"]
        ) as fa:
            if kwargs['full_header']:
                _fasta_with_full_header_assets(fa)
            else:
                _fasta_without_full_header_assets(fa)


@pytest.mark.parametrize(
    argnames="kwargs",
    argvalues=(
            {
                "newline": test_case[0],
                "fasta_seq": test_case[1]
            } for test_case in itertools.product(
        VALID_NEWLINE,
        [test_oneline_fasta_seq, test_fasta_seq]
    )
    )
)
def test_fai(kwargs) -> None:
    try:
        from pysam import faidx
    except (ImportError, ValueError): # ValueError raised by pypy
        return
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_filename = os.path.join(tmpdir, "tmp.fa")
        with get_writer(fasta_filename, newline=kwargs["newline"]) as fh:
            fh.write(kwargs["fasta_seq"])
        faidx(fasta_filename)
        try:
            create_fai_from_fasta(fasta_filename, fasta_filename + ".tetgs.fai")
        except FastaIndexNotWritableError:
            return
        assert open(fasta_filename + ".fai").read() == open(fasta_filename + ".tetgs.fai").read()


def test_subset_to_file():
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_filename = os.path.join(tmpdir, "tmp.fa")
        with get_writer(fasta_filename) as fh:
            fh.write(test_fasta_seq)
        with FastaViewFactory(
                filename=fasta_filename,
                full_header=False,
                read_into_memory=True
        ) as fa:
            subset_filename = fasta_filename + ".subset.fa"
            querys = (('chr3', 2, 15), ('chr3', 2), ('chr3',))
            with pytest.raises(DuplicatedChromosomeNameError):
                fa.subset_to_file(subset_filename, querys)
            fa.subset_to_file(subset_filename, querys, ["1", "2", "3"])
            # TODO: Not finished
