import os
import random
import string
import tempfile

import pytest

from labw_utils.commonutils.lwio import get_reader, get_writer, TqdmReaderProxy, TqdmLineReaderProxy
from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_reader, get_tqdm_line_reader
from labw_utils.commonutils.stdlib_helper import shutil_helper
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf

CONTENT_SIZE = 1024


def assess_binary_archive_io(filename: str):
    available_chars = string.printable

    contents = bytes(
        "".join((random.choices(available_chars, k=CONTENT_SIZE))),
        encoding="UTF-8"
    )
    len_contents = len(contents)

    shutil_helper.rm_rf(filename)
    with get_writer(filename, is_binary=True, compression_level=9) as writer:
        writer.write(contents)
    with get_reader(filename, is_binary=True) as reader:
        assert reader.read(len_contents) == contents
    with get_tqdm_reader(filename, is_binary=True) as reader:
        reader: TqdmReaderProxy
        assert reader.read(len_contents) == contents
        assert reader._tqdm._total == len(contents)
    assert filename.endswith("txt") == (os.path.getsize(filename) == len_contents)


def assess_text_archive_io(filename: str):
    available_chars = string.ascii_letters + "\n"

    contents = "".join(random.choices(available_chars, k=CONTENT_SIZE))
    len_contents = len(contents)

    contents_list = contents.split('\n')
    shutil_helper.rm_rf(filename)
    with get_writer(filename, is_binary=False, newline='\n', compression_level=9) as writer:
        writer.write(contents)
    with get_reader(filename, is_binary=False, newline='\n') as reader:
        assert reader.read(len_contents) == contents
    with get_tqdm_reader(filename, newline='\n') as reader:
        reader: TqdmReaderProxy
        assert reader.read(len_contents) == contents
        assert reader._tqdm._total == len(contents)
    with get_tqdm_line_reader(filename, newline='\n') as reader:
        reader: TqdmLineReaderProxy
        i = 0
        assert reader._tqdm._total == len(contents_list)
        for line in reader:
            assert contents_list[i] == line
            i += 1
            assert reader._tqdm._n == i


extensions = (
    "txt", "xz", "bz2", "lzma", "gz"
)


@pytest.mark.parametrize(
    argnames="ext",
    argvalues=extensions,
    ids=["test_" + ext for ext in extensions]
)
def test_ext(ext: str):
    with tempfile.TemporaryDirectory() as tempdir:
        filename = os.path.join(tempdir, f"1.{ext}")
        assess_text_archive_io(filename)
        assess_binary_archive_io(filename)
        rm_rf(filename)
