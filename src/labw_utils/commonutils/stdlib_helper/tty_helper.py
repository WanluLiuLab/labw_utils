"""
``labw_utils.stdlib_helper.tty_helper`` -- ANSI colors

.. note ::
    Please make sure that the file descriptor you wish to apply is a Tele-TypeWriter (TTY)
    before applying ANSI colors! Convenient checking by :py:func:`os.isatty`.

.. versionadded:: 1.0.2
"""

__all__ = ("STANDARD_ANSI_OPERATORS", "get_ansi_rgb")

STANDARD_ANSI_OPERATORS = {
    "RESET": "\033[0m",
    "BOLD": "\033[1m",
    "REVERSE": "\033[7m",
    "UNDERLINE": "\033[4m",
    "NO_BOLD": "\033[22m",
    "NO_UNDERLINE": "\033[24m",
    "NO_REVERSE": "\033[27m",
    "FG_RED": "\033[31m",
    "FG_GREEN": "\033[32m",
    "FG_YELLOW": "\033[33m",
    "FG_BLUE": "\033[34m",
    "FG_PURPLE": "\033[35m",
    "FG_CRAYON": "\033[36m",
    "FG_WHITE": "\033[37m",
    "FG_DEFAULT": "\033[39m",
    "BG_RED": "\033[41m",
    "BG_GREEN": "\033[42m",
    "BG_YELLOW": "\033[43m",
    "BG_BLUE": "\033[44m",
    "BG_PURPLE": "\033[45m",
    "BG_CRAYON": "\033[46m",
    "BG_WHITE": "\033[47m",
    "BG_DEFAULT": "\033[49m",
}
"""
A subset of standard ANSI operations

.. versionadded:: 1.0.2
"""


def get_ansi_rgb(r: int, g: int, b: int, is_fg: bool = True) -> str:
    """
    Create ANSI color in RGB mode.

    :param r: Red color, 0 to 255.
    :param g: Green color, 0 to 255.
    :param b: Blue color, 0 to 255.
    :param is_fg: Whether it is frontend or backend.
    :return: Generated ANSI control sequence.

    .. versionadded:: 1.0.2
    """
    if is_fg:
        return f"\x1b[38;2;{r};{g};{b}m"
    else:
        return f"\x1b[48;2;{r};{g};{b}m"
