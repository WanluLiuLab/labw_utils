import os
import sqlite3
import tempfile
import time

import pandas as pd
import pytest

from labw_utils.commonutils.appender import load_table_appender_class, list_table_appender, TableAppenderConfig, \
    BaseTableAppender
from labw_utils.commonutils.stdlib_helper.shutil_helper import wc_l

AVAILABLE_TABLE_APPENDERS = [v_doc[0] for v_doc in list_table_appender()]


def validate_lines(appender: BaseTableAppender, required_number_of_lines: int) -> None:
    actual_number_of_lines = 0
    if appender.real_filename == "":
        return
    elif appender.real_filename.endswith("hdf5"):
        actual_number_of_lines = pd.read_hdf(appender.real_filename, key="df").shape[0]
    elif "tsv" in appender.real_filename:
        actual_number_of_lines = wc_l(appender.real_filename) - 1
    elif appender.real_filename.endswith("parquet"):
        actual_number_of_lines = pd.read_parquet(appender.real_filename).shape[0]
    elif appender.real_filename.endswith("sqlite3"):
        with sqlite3.connect(appender.real_filename) as con:
            actual_number_of_lines = pd.read_sql_query("SELECT * FROM db", con=con).shape[0]
    if actual_number_of_lines != required_number_of_lines:
        raise AssertionError(
            f"{appender.real_filename}, "
            f"Required: {required_number_of_lines} "
            f"Actual: {actual_number_of_lines}"
        )


@pytest.mark.parametrize(
    argnames="name",
    argvalues=AVAILABLE_TABLE_APPENDERS,
    ids=["test_" + appender_name for appender_name in AVAILABLE_TABLE_APPENDERS]
)
def test_appender(name: str):
    with tempfile.TemporaryDirectory() as tmpdir:
        for lines_to_append in (1, 4, 5, 6, 9, 10, 11, 1024):
            with load_table_appender_class(name)(
                    filename=os.path.join(tmpdir, "test"),
                    header=("INDEX", "TIME", "3"),
                    tac=TableAppenderConfig(buffer_size=5)
            ) as appender:
                for i in range(lines_to_append):
                    appender.append((i, time.asctime(), "3"))
            validate_lines(appender, lines_to_append)
            if appender.real_filename != '':
                os.remove(appender.real_filename)
