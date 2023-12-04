"""l
ist_appenders -- List all available table appenders.

.. versionadded:: 1.0.2
"""
from typing import List

from labw_utils.commonutils.appender import list_table_appender


def main(_: List[str]):
    for appender in list_table_appender():
        print(": ".join(appender))
    return 0
