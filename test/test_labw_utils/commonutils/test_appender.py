import os

import pytest

from labw_utils.commonutils.appender import load_table_appender_class, AVAILABLE_TABLE_APPENDERS
from labw_utils.commonutils.appender.typing import TableAppenderConfig, BaseTableAppender


def assert_appender(
        appender: BaseTableAppender,
        lines_to_append: int
):
    for i in range(lines_to_append):
        appender.append([i, "2", "3"])
    appender.close()
    appender.validate_lines(lines_to_append)


@pytest.mark.parametrize(
    argnames="name",
    argvalues=AVAILABLE_TABLE_APPENDERS,
    ids=["test_" + appender_name for appender_name in AVAILABLE_TABLE_APPENDERS]
)
def test_appender(name: str):
    for lines_to_append in (1, 4, 5, 6, 9, 10, 11, 1024):
        appender = load_table_appender_class(name)(
            filename="test",
            header=["1", "2", "3"],
            tac=TableAppenderConfig(buffer_size=5)
        )
        assert_appender(appender, lines_to_append)
        if appender.real_filename != '':
            os.remove(appender.real_filename)
