import os
import random
import string

import pytest

import conftest
from labw_utils.commonutils.io.safe_io import get_writer, get_reader
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_reader, get_tqdm_line_reader, TqdmReader
from labw_utils.commonutils.stdlib_helper import shutil_helper


@pytest.fixture(scope="module", autouse=True)
def initialize_module(initialize_session) -> conftest.ModuleTestInfo:
    """
    This function sets up a directory for testing
    """
    session_test_info = initialize_session
    module_test_info = conftest.ModuleTestInfo(session_test_info.base_test_dir, __name__)
    yield module_test_info
    module_test_info.teardown()


CONTENT_SIZE = 1024


def assess_binary_archive_io(filename: str):
    available_chars = string.printable

    contents = bytes("".join((random.choice(available_chars) for _ in range(CONTENT_SIZE))), encoding="UTF-8")
    len_contents = len(contents)

    shutil_helper.rm_rf(filename)
    with get_writer(filename, is_binary=True) as writer:
        writer.write(contents)
    with get_reader(filename, is_binary=True) as reader:
        assert reader.read(len_contents) == contents
    with get_tqdm_reader(filename, is_binary=True) as reader:
        reader: TqdmReader
        assert reader.read(len_contents) == contents
        assert reader._tqdm.total == len(contents)
    assert filename.endswith("txt") == (os.path.getsize(filename) == len_contents)


def assess_text_archive_io(filename: str):
    # FIXME: Bugs here
    available_chars = string.printable

    contents = "".join((random.choice(available_chars) for _ in range(CONTENT_SIZE))).replace('\r', '')
    len_contents = len(contents)

    contents_list = contents.split('\n')
    shutil_helper.rm_rf(filename)
    with get_writer(filename, encoding="UTF-8", newline='\n') as writer:
        writer.write(contents)
    with get_reader(filename, encoding="UTF-8", newline='\n') as reader:
        assert reader.read(len_contents) == contents
    with get_tqdm_reader(filename, encoding="UTF-8", newline='\n') as reader:
        reader: TqdmReader
        assert reader.read(len_contents) == contents
        assert reader._tqdm.total == len(contents)
    with get_tqdm_line_reader(filename, encoding="UTF-8", newline='\n') as reader:
        reader: TqdmReader
        i = 0
        assert reader._tqdm.total == len(contents_list)
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
def test_ext(initialize_module, ext: str):
    module_test_info = initialize_module
    filename = os.path.join(module_test_info.path, f"1.{ext}")
    assess_text_archive_io(filename)
    assess_binary_archive_io(filename)
    shutil_helper.rm_rf(filename)
