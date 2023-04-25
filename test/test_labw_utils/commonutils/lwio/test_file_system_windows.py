import os
import pathlib

import pytest

from labw_utils.commonutils.lwio import file_system

if os.name != 'nt':
    pytest.skip("System is not NT", allow_module_level=True)


def test_get_abspath():
    assert file_system.get_abspath('') == ''
    home = str(pathlib.Path.home())
    assert file_system.get_abspath('~') == home
    assert file_system.get_abspath('~\\..\\') == os.path.dirname(home)
    assert file_system.get_abspath('~\\~') == home + '\\~'


def test_file_exists():
    pass  # TODO


def test_directory_exists():
    pass  # TODO


def test_is_soft_link():
    pass  # TODO
