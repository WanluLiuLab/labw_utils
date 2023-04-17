from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List

blh = get_logger("b")
ulh = get_logger()


def main(_: List[str]):
    blh.trace("BLH TRACE")
    blh.debug("BLH DEBUG")
    blh.info("BLH INFO")
    blh.warning("BLH WARNING")

    ulh.trace("ULH TRACE")
    ulh.debug("ULH DEBUG")
    ulh.info("ULH INFO")
    ulh.warning("ULH WARNING")
