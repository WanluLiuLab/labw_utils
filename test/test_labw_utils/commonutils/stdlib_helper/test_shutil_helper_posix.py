import os

import pytest

import conftest
from labw_utils.commonutils.lwio.file_system import get_abspath
from labw_utils.commonutils.stdlib_helper import shutil_helper

if os.name != "posix":
    pytest.skip("System is NOT POSIX", allow_module_level=True)


@pytest.fixture(scope="module", autouse=True)
def initialize_module(initialize_session) -> conftest.ModuleTestInfo:
    """
    This function sets up a directory for testing
    """
    session_test_info = initialize_session
    module_test_info = conftest.ModuleTestInfo(session_test_info.base_test_dir, __name__)
    yield module_test_info
    module_test_info.teardown()


def test_wc_c():
    assert shutil_helper.wc_c("/dev/null") == 0
    assert shutil_helper.wc_c("/dev/zero") == 0
    assert shutil_helper.wc_c("/dev/stdin") == 0


def test_readlink_f(initialize_module):
    assert shutil_helper.readlink_f("") == ""
    test_path = initialize_module.path
    if os.name == "posix":
        shutil_helper.rm_rf(test_path)
        shutil_helper.touch(os.path.join(test_path, "aa"))
        os.symlink(os.path.join(test_path, "aa"), os.path.join(test_path, "ab"))
        os.symlink(os.path.join(test_path, "ab"), os.path.join(test_path, "ac"))
        assert shutil_helper.readlink_f(os.path.join(test_path, "ab")) == get_abspath(os.path.join(test_path, "aa"))
        assert shutil_helper.readlink_f(get_abspath(os.path.join(test_path, "ab"))) == get_abspath(
            os.path.join(test_path, "aa")
        )
        assert shutil_helper.readlink_f(os.path.join(test_path, "ac")) == get_abspath(os.path.join(test_path, "aa"))
        assert shutil_helper.readlink_f(get_abspath(os.path.join(test_path, "ac"))) == get_abspath(
            os.path.join(test_path, "aa")
        )
        shutil_helper.rm_rf(test_path)
