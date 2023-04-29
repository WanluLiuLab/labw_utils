from pyspark.sql import SparkSession

from labw_utils.commonutils.libfrontend import setup_basic_logger
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)
setup_basic_logger()
spark = SparkSession.builder.getOrCreate()
_lh.info("Using pyspark version %s", spark.version)


