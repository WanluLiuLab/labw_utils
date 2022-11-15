import logging

import pytest

from labw_utils.commonutils.stdlib_helper import logger_helper


def test_logger():
    lh = logging.getLogger()
    with pytest.raises(ValueError):
        logger_helper.set_level('NON_EXISTENT_LEVEL')
    lh.debug("DEBUG")
    lh.info("INFO")
    lh.warning("WARNING")
    lh.error("ERROR")
    logger_helper.set_level(logging.DEBUG)
