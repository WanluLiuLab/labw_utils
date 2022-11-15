import os

import pytest

from labw_utils.commonutils.stdlib_helper import shutil_helper

if os.name != 'nt':
    pytest.skip("System is NOT NT", allow_module_level=True)


def test_wc_c(initialize_module):
    assert shutil_helper.wc_c('nul') == 0
