"""
``labw_utils.stdlib_helper.logger_helper`` -- Additional logging facilities.

It performs the following:

- Addition of a ``TRACE`` level, which is ``8`` and less than ``DEBUG``.
- An easy-to-use logger setup facility.

.. versionadded:: 1.0.2
"""

import logging
import sys
from logging import DEBUG, WARNING, ERROR, FATAL, INFO

from labw_utils.typing_importer import Optional, Union

__all__ = ("DEBUG", "WARNING", "ERROR", "FATAL", "INFO", "TRACE", "get_logger", "get_formatter")

TRACE = 8
"""
New logging level

.. versionadded:: 1.0.2
"""


def trace(self, msg, *args, **kwargs):
    """
    Log 'msg % args' with severity 'TRACE'.
    """
    if self.isEnabledFor(TRACE):
        self._log(TRACE, msg, args, **kwargs)


logging.addLevelName(TRACE, "TRACE")
logging.Logger.trace = trace  # type: ignore
logging.trace = trace  # type: ignore


def get_formatter(level: Union[int, str]) -> logging.Formatter:
    """
    Create standard formatter.

    :param level: Cutoff level of the logger.
    :return: Generated formatter. Will be:
      - ``%(asctime)s\t[%(levelname)s] %(message)s`` if level is larger than debug.
      - ``%(asctime)s %(name)s:%(lineno)d::%(funcName)s\t[%(levelname)s]\t%(message)s`` otherwise.
    :raises ValueError: If level does not exist.

    .. versionadded:: 1.0.2
    """
    if isinstance(level, str):
        level = logging.getLevelName(level)
    if isinstance(level, str):
        raise ValueError(f"{level} not exist!")
    if level > logging.DEBUG:
        log_format = "%(asctime)s\t[%(levelname)s] %(name)s\t%(message)s"
    else:
        log_format = "%(asctime)s %(name)s:%(lineno)d::%(funcName)s\t[%(levelname)s]\t%(message)s"
    return logging.Formatter(log_format)


def get_logger(
    name: Optional[str] = None,
    level: Union[str, int] = TRACE,
    log_to_stderr: bool = False,
    log_stderr_level: Union[str, int] = INFO,
    log_stderr_formatter: Optional[logging.Formatter] = None,
    log_file_name: Optional[str] = None,
    log_file_level: Union[str, int] = TRACE,
    log_file_formatter: Optional[logging.Formatter] = None,
) -> logging.Logger:
    """
    A Simple :py:func:`logging.getLogger()` wrapper.

    :param name: Name of the logger. Recommended to be module name (i.e., ``__name__`` global in module).
    :param level: Cutoff level of the logger.
    :param log_to_stderr: Whether to log to standard error.
    :param log_stderr_level: Standard error logging level.
    :param log_stderr_formatter: Standard error logging formatter. If not set, would automatically generate
      one from ``log_stderr_level``.
    :param log_file_name: Filename of the log. Set to :py:obj:`None` to avoid creating files.
    :param log_file_level: File logging level.
    :param log_file_formatter: File logging formatter. If not set, would automatically generate
      one from ``log_stderr_level``.
    :return: The logger handler.
    :raises ValueError: If level does not exist.

    .. versionadded:: 1.0.2
    """
    if name is None:
        return logging.getLogger()
    else:
        logger = logging.getLogger(name)
        logger.setLevel(level)
        if log_file_name is not None:
            file_handler = logging.FileHandler(filename=log_file_name)
            file_handler.setLevel(log_file_level)
            if log_file_formatter is not None:
                file_handler.setFormatter(log_file_formatter)
            else:
                file_handler.setFormatter(get_formatter(log_file_level))
            logger.addHandler(file_handler)
        if log_to_stderr:
            serr_handler = logging.StreamHandler(stream=sys.stderr)
            serr_handler.setLevel(log_stderr_level)
            if log_stderr_formatter is not None:
                serr_handler.setFormatter(log_stderr_formatter)
            else:
                serr_handler.setFormatter(get_formatter(log_stderr_level))
            logger.addHandler(serr_handler)
        return logger
