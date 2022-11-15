import os
import sys
import tempfile

import pytest

import conftest
from labw_utils.commonutils.io.file_system import file_exists
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


def test_mkdir_p_and_rm_f(initialize_module):
    aa = os.path.join(initialize_module.path, "aa")
    shutil_helper.touch(aa)
    with pytest.raises(IsADirectoryError):
        shutil_helper.touch(initialize_module.path)
    assert os.path.isdir(initialize_module.path)
    assert os.path.isfile(aa)
    shutil_helper.rm_rf(aa)
    shutil_helper.touch(aa)
    shutil_helper.rm_rf(initialize_module.path)
    assert not os.path.isdir(initialize_module.path)
    assert not os.path.isfile(aa)


def test_wc_c(initialize_module):
    aa = os.path.join(initialize_module.path, "aa")
    shutil_helper.touch(aa)
    assert shutil_helper.wc_c(aa) == 0
    shutil_helper.rm_rf(aa)
    self_exe = sys.executable
    if file_exists(self_exe):
        assert shutil_helper.wc_c(self_exe) == os.path.getsize(self_exe)


@pytest.mark.parametrize(
    argnames="kwargs",
    argvalues=(
            {"text": b"", "answer": 0},
            {"text": b"\n", "answer": 1},
            {"text": b"A", "answer": 1},
            {"text": b"A\n", "answer": 1},
    )
)
def test_wc_l(kwargs):
    with tempfile.TemporaryDirectory() as tmpdir:
        target_filepath = os.path.join(tmpdir, "tmp")
        with open(target_filepath, "wb") as writer:
            writer.write(kwargs["text"])
        assert shutil_helper.wc_l(target_filepath) == kwargs["answer"]
