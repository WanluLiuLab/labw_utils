"""
logger_helper.py -- System-Wide Logger.

Features
--------

This logger defines a "trace" level (TRACE = 8) and a logging decorator.
It can also read environment variable named 'LOG_LEVEL' and fallback to DEBUG (10) by default.

Usage
-----

See :py:func:`chronolog` and :py:func:`trace`.

There are also :py:func:`set_level` and :py:func:`get_logger`,
which is only snake-case wrappers for those contents inside :py:mod:`logging` standard module.
"""

import logging
import os
import sys
from logging import DEBUG, WARNING, ERROR, FATAL
from typing import Optional, Union

__all__ = (
    'TRACE',
    'DEBUG',
    'WARNING',
    'ERROR',
    'FATAL',
    'chronolog',
    'get_logger'
)

_SB = os.environ.get('SPHINX_BUILD')

# The following contents adds a new log level called trace.
TRACE = 8


def trace(self, msg, *args, **kwargs):
    """
    Log 'msg % args' with severity 'TRACE'.
    """
    if self.isEnabledFor(TRACE):
        self._log(TRACE, msg, args, **kwargs)


logging.addLevelName(TRACE, "TRACE")
logging.Logger.trace = trace
logging.trace = trace

_lh = logging.getLogger(__name__)


def get_formatter(level: Union[int, str]) -> logging.Formatter:
    if isinstance(level, str):
        level = logging.getLevelName(level)
    if isinstance(level, str):
        raise ValueError(f"{level} not exist!")
    if level > logging.DEBUG:
        log_format = '%(asctime)s\t[%(levelname)s] %(message)s'
    else:
        log_format = '%(asctime)s %(name)s:%(lineno)d::%(funcName)s\t[%(levelname)s]\t%(message)s'
    return logging.Formatter(log_format)


def set_level(**_):
    """
    Deprecated dumb function
    """
    pass


def chronolog(display_time: bool = False, log_error: bool = False):
    """
    The logging decorator, will inject a logger variable named _lh to the code.
    From <https://stackoverflow.com/questions/17862185/how-to-inject-variable-into-scope-with-a-decorator>

    .. note::
        The :py:func:`error` (or :py:func:`exception`, :py:func:`critical`, :py:func:`fatal`
        functions DO NOT exit the program! You have to exit the program by yourself!

    .. warning::
        Call this function, do NOT call functions inside this function!

    :param display_time: Whether to display calling time, arguments and return value in [TRACE] log level.
    :param log_error: Whether add error captured
    """

    def msg_decorator(f):
        if _SB == '1':
            return f  # To make Sphinx get the right result.

        def inner_dec(*args, **kwargs):
            """
            Decorator which performs the logging and do the work.

            :param args: Unnamed arguments of the decorated function call.
            :param kwargs: Named arguments of the decorated function call.
            :return: The return value of the decorated function call.
            :raise: The return value of the decorated function call.
            """
            try:
                _ = f.__globals__
            except AttributeError:
                return f(*args, **kwargs)
            lh = logging.getLogger(f.__module__)
            if display_time:
                args_repr = [repr(a) for a in args]
                kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
                signature = ", ".join(args_repr + kwargs_repr)
                # FIXME: Function name incorrect!
                lh.trace(f"{f.__name__}({signature})", )
            res = None
            try:
                res = f(*args, **kwargs)
            except Exception as e:
                if log_error:
                    lh.exception(f"Exception inside func: {e}", stack_info=True, exc_info=True)
                raise e
            finally:
                if display_time:
                    lh.trace(f"{f.__name__} -> {res}", )
            return res

        return inner_dec

    return msg_decorator


def get_logger(
        name: Optional[str] = None,
        level: Optional[Union[str, int]] = TRACE,
        log_to_stderr: bool = False,
        log_stderr_level: Optional[Union[str, int]] = logging.INFO,
        log_stderr_formatter: Optional[logging.Formatter] = None,
        log_file_name: Optional[str] = None,
        log_file_level: Optional[Union[str, int]] = TRACE,
        log_file_formatter: Optional[logging.Formatter] = None
):
    """
    A Simple logging.getLogger() wrapper.
    """
    if name is None:
        return logging.getLogger()
    else:
        logger = logging.getLogger(name)
        logger.setLevel(level)
        if log_file_name is not None:
            file_handler = logging.FileHandler(
                filename=log_file_name
            )
            file_handler.setLevel(log_file_level)
            if log_file_formatter is not None:
                file_handler.setFormatter(log_file_formatter)
            else:
                file_handler.setFormatter(get_formatter(log_file_level))
            logger.addHandler(
                file_handler
            )
        if log_to_stderr:
            serr_handler = logging.StreamHandler(
                stream=sys.stderr
            )
            serr_handler.setLevel(log_stderr_level)
            if log_stderr_formatter is not None:
                serr_handler.setFormatter(log_stderr_formatter)
            else:
                serr_handler.setFormatter(get_formatter(log_stderr_level))
            logger.addHandler(
                serr_handler
            )
        return logger
