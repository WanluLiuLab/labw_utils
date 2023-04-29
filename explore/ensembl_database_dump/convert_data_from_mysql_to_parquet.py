"""
convert_data_from_mysql_to_parquet.py -- Synchronise data from Ensembl to local database and Apache Parquet file.
"""

import os
import pandas as pd
import pyarrow as pa
import sqlalchemy
from labw_utils.commonutils.libfrontend import setup_basic_logger
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from pyarrow import parquet as pq

# Constants

DB_NAME = "homo_sapiens_core_109_38"

# Services
setup_basic_logger()

_lh = get_logger(__name__)


def read_table_to_parquet(table_name: str):
    _flh = get_logger(__name__)
    _engine = sqlalchemy.create_engine(
        f"mariadb+pymysql://ensdb:ensdb@localhost:13306/{DB_NAME}?charset=utf8mb4"
    )
    with _engine.connect() as con:
        _flh.info("Reading %s", table_name)
        table_contents = pd.read_sql_table(
            table_name=table_name,
            con=con,
            dtype_backend="pyarrow"
        )

    _flh.info("Writing %s", table_name)
    pq.write_table(
        pa.Table.from_pandas(df=table_contents),
        os.path.join("converted_parquet", table_name + ".parquet"),
        flavor='spark'
    )
    _flh.info("Table %s FIN", table_name)


def read_db_to_parquet():
    _lh.info("Getting table names...")
    os.makedirs("converted_parquet", exist_ok=True)
    engine = sqlalchemy.create_engine(
        f"mariadb+pymysql://ensdb:ensdb@localhost:13306/{DB_NAME}?charset=utf8mb4"
    )
    engine_inspector = sqlalchemy.inspect(engine)
    table_names = engine_inspector.get_table_names()
    _lh.info("Got %d tables", len(table_names))

    # Single-Threaded due to memory usage on large datasets.
    for table_name in table_names:
        read_table_to_parquet(table_name)


if __name__ == '__main__':
    setup_basic_logger()
    read_db_to_parquet()
