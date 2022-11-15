import os
from enum import Enum


class AnsiColorEnum(Enum):
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    PURPLE = "\033[35m"
    CRAYON = "\033[36m"
    CLEAR = "\033[0m"


class DumbColorEnum(Enum):
    RED = ""
    GREEN = ""
    YELLOW = ""
    BLUE = ""
    PURPLE = ""
    CRAYON = ""
    CLEAR = ""


def get_color(fd: int):
    """
    Get ANSI color dictionary for current file descriptor.

    Python can test whether the output is a tty. Other method have to employ ncurses.

    :param fd: File descriptor.
    :return: A color enum with a format like ``{RED = "\\033[31m"}``.
    """
    if os.isatty(fd):
        retd = AnsiColorEnum
    else:
        retd = DumbColorEnum
    return retd
