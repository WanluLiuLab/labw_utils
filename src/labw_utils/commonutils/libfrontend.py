"""
``labw_utils.commonutils.libfrontend`` -- Helpers to setup commandline frontend that have multiple subcommands.

.. versionadded:: 1.0.2

Consumed environment variables:

.. envvar:: LOG_LEVEL

can either be string or integer.

.. envvar:: LOG_FILE_NAME


"""

from __future__ import annotations

__all__ = (
    "setup_frontend",
    "get_subcommands",
    "get_argparser_from_subcommand",
    "get_main_func_from_subcommand",
    "get_doc_from_subcommand",
    "setup_basic_logger",
    "add_file_handler_to_root_logger_handler",
)

import argparse
import importlib
import inspect
import logging
import os
import pkgutil
import sys

from labw_utils import UnmetDependenciesError
from labw_utils.commonutils.stdlib_helper import logger_helper
from labw_utils.stdlib.cpy310.pkgutil import resolve_name
from labw_utils.typing_importer import Iterable
from labw_utils.typing_importer import Optional, Callable, List

_lh: Optional[logging.Logger] = None


def setup_basic_logger():
    """
    Clear current logging configuration and setup basic command-line logging support.
    Would respect :envvar:`LOG_LEVEL` environment variable.

    .. versionadded:: 1.0.2
    """
    _stream_handler = logging.StreamHandler()
    _stream_handler.setLevel(os.environ.get("LOG_LEVEL", "INFO"))
    _stream_handler.setFormatter(logger_helper.get_formatter(_stream_handler.level))

    if sys.version_info < (3, 8):
        logging.basicConfig(handlers=[_stream_handler], level=logger_helper.TRACE)
    else:
        logging.basicConfig(handlers=[_stream_handler], force=True, level=logger_helper.TRACE)


_NONE_DOC = "NONE DOC"


def get_subcommands(package_main_name: str, verbose: bool = False) -> Iterable[str]:
    """
    Get valid name of subcommands.
    Would automatically skip import failures.

    .. versionadded:: 1.0.2
    """
    for spec in pkgutil.iter_modules(resolve_name(package_main_name).__spec__.submodule_search_locations):
        subcommand_name = spec.name
        if subcommand_name.startswith("_"):
            continue
        try:
            _ = resolve_name(f"{package_main_name}.{subcommand_name}")
            _ = get_main_func_from_subcommand(package_main_name, subcommand_name)
        except (UnmetDependenciesError, ImportError, AttributeError):
            if verbose and _lh is not None:
                _lh.warning("Subcommand %s have unmet dependencies!", subcommand_name)
            continue
        yield subcommand_name


def get_doc_from_subcommand(package_main_name: str, subcommand_name: str) -> Optional[str]:
    """
    Return documentation of that module

    .. versionadded:: 1.0.2
    """
    importlib.import_module(f"{package_main_name}.{subcommand_name}")
    i = resolve_name(f"{package_main_name}.{subcommand_name}")
    if hasattr(i, "__doc__"):
        return i.__doc__
    else:
        return None


def get_main_func_from_subcommand(package_main_name: str, subcommand_name: str) -> Optional[Callable[[list[str]], int]]:
    """
    Return a subcommands' :py:func:`main` function.

    .. versionadded:: 1.0.2
    """
    importlib.import_module(f"{package_main_name}.{subcommand_name}")
    i = resolve_name(f"{package_main_name}.{subcommand_name}")
    if hasattr(i, "main") and inspect.isfunction(getattr(i, "main")):
        return i.main
    else:
        return None


def get_argparser_from_subcommand(package_main_name: str, subcommand_name: str) -> Optional[argparse.ArgumentParser]:
    """
    Return result of a subcommands' :py:func:`create_parser` function.

    .. versionadded:: 1.0.2
    """
    importlib.import_module(f"{package_main_name}.{subcommand_name}")
    i = resolve_name(f"{package_main_name}.{subcommand_name}")
    if hasattr(i, "create_parser") and inspect.isfunction(getattr(i, "create_parser")):
        return i.create_parser()
    else:
        return None


def lscmd(package_main_name: str, valid_subcommand_names: Iterable[str]):
    """
    ``lscmd`` frontend

    .. versionadded:: 1.0.2
    """
    name_doc_dict = {}
    if _lh is not None:
        _lh.info("Listing modules...")
    for item in valid_subcommand_names:
        doc = get_doc_from_subcommand(package_main_name, item)
        if doc is None:
            doc = _NONE_DOC
        else:
            doc_splitlines = doc.splitlines()
            if not doc_splitlines:
                doc = _NONE_DOC
            else:
                while len(doc_splitlines) > 0:
                    potential_doc = doc_splitlines[0].strip()
                    if potential_doc == "":
                        doc_splitlines.pop(0)
                    else:
                        doc = potential_doc
                        break
                else:
                    doc = _NONE_DOC
        if doc.find("--") != -1:
            doc = doc.split("--")[1].strip()
        name_doc_dict[item] = doc
    sys.stderr.flush()
    sys.stdout.flush()
    for item, doc in name_doc_dict.items():  # To prevent logger from polluting outout
        print(f"{item} -- {doc}")
    sys.exit(0)


class _ParsedArgs:
    input_subcommand_name: str = ""
    have_help: bool = False
    have_version: bool = False
    parsed_args: List[str] = []


def _parse_args(args: List[str]) -> _ParsedArgs:
    parsed_args = _ParsedArgs()
    i = 0
    while i < len(args):
        name = args[i]
        if name in ("--help", "-h"):
            parsed_args.have_help = True
        elif name in ("--version", "-v"):
            parsed_args.have_version = True
        elif not name.startswith("-") and parsed_args.input_subcommand_name == "":
            parsed_args.input_subcommand_name = name
            args.pop(i)
        i += 1
    parsed_args.parsed_args = args
    return parsed_args


def _format_help_info(prefix: str, package_main_name: str) -> str:
    return f"""
This is frontend of `{package_main_name.split('.')[0].strip()}` provided by `{__name__}`.

SYNOPSYS: {prefix} [[SUBCOMMAND] [ARGS_OF SUBCOMMAND] ...] [-h|--help] [-v|--version]

If a valid [SUBCOMMAND] is present, will execute [SUBCOMMAND] with all other arguments

If no valid [SUBCOMMAND] is present, will fail to errors.

If no [SUBCOMMAND] is present, will consider options like:
    [-h|--help] show this help message and exit
    [-v|--version] show package version and other information

ENVIRONMENT VARIABLES:

    LOG_FILE_NAME: The sink of all loggers.
    LOG_LEVEL: Targeted frontend log level. May be DEBUG INFO WARN ERROR FATAL.

Use `lscmd` as subcommand with no options to see available subcommands.
"""


def add_file_handler_to_root_logger_handler(default_log_filename: str):
    """
    Register a file handler.

    Wouls respect :envvar:`LOG_FILE_NAME` environment variable.

    .. versionadded:: 1.0.2
    """
    log_filename = os.environ.get("LOG_FILE_NAME", default_log_filename)
    file_handler = logging.FileHandler(filename=log_filename)
    file_handler.setLevel(logger_helper.TRACE)
    file_handler.setFormatter(logger_helper.get_formatter(logger_helper.TRACE))
    logging.root.addHandler(file_handler)


def setup_frontend(
    package_main_name: str,
    one_line_description: str,
    version: str,
    help_info: Optional[str] = None,
    subcommand_help: str = "Use 'lscmd' to list all valid subcommands.",
    use_root_logger: bool = True,
    default_log_filename: Optional[str] = None,
):
    """
    Setup the frontend.

    :param package_main_name: TODO
    :param one_line_description: TODO
    :param version: TODO
    :param help_info: TODO
    :param subcommand_help: TODO
    :param use_root_logger: TODO
    :param default_log_filename: TODO

    .. versionadded:: 1.0.2
    """
    global _lh
    if default_log_filename is None:
        pkg_full_name = package_main_name.split(".")
        if pkg_full_name[-1] in {"main", "_main", "__main__"}:
            pkg_full_name.pop()
        default_log_filename = ".".join(pkg_full_name) + ".log"
    setup_basic_logger()
    _lh = logger_helper.get_logger(__name__)
    _lh.info(f"{one_line_description} ver. {version}")
    argv = sys.argv
    _lh.info(f'Called by: {" ".join(argv)}')
    parsed_args = _parse_args(argv[1:])
    if use_root_logger:
        add_file_handler_to_root_logger_handler(default_log_filename)
    if parsed_args.input_subcommand_name == "lscmd":
        get_subcommands_verbose = True
    else:
        get_subcommands_verbose = False
    valid_subcommand_names = get_subcommands(package_main_name, get_subcommands_verbose)
    if parsed_args.input_subcommand_name == "lscmd":
        lscmd(package_main_name, valid_subcommand_names)
    elif parsed_args.input_subcommand_name == "":
        if parsed_args.have_help:
            if help_info is None:
                help_info = _format_help_info(argv[0], package_main_name)
            print(help_info)
            sys.exit(0)
        elif parsed_args.have_version:
            print(version)
            sys.exit(0)
        else:
            _lh.exception(f"Subcommand name not set! {subcommand_help}")
            sys.exit(1)
    else:
        main_fnc = get_main_func_from_subcommand(
            package_main_name=package_main_name, subcommand_name=parsed_args.input_subcommand_name
        )
        if main_fnc is not None:
            sys.exit(main_fnc(parsed_args.parsed_args))
        else:
            _lh.exception(f"Subcommand '{parsed_args.input_subcommand_name}' not found! {subcommand_help}")
            sys.exit(1)
