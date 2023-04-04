__all__ = (
    'setup_frontend',
    "get_subcommands",
    "get_argparser_from_subcommand",
    "get_main_func_from_subcommand"
)

import argparse
import importlib
import inspect
import logging
import os
import pkgutil
import sys
from typing import List, Iterable, Callable, Optional

from labw_utils.commonutils.stdlib_helper import logger_helper
from labw_utils.stdlib.cpy310.pkgutil import resolve_name

stream_handler = logging.StreamHandler()
stream_handler.setLevel(os.environ.get('LOG_LEVEL', 'INFO'))
stream_handler.setFormatter(logger_helper.get_formatter(stream_handler.level))

logging.basicConfig(
    handlers=[
        stream_handler
    ],
    force=True,
    level=logger_helper.TRACE
)
_lh = logger_helper.get_logger(__name__)

NONE_DOC = "NONE DOC"


def get_subcommands(package_main_name: str) -> Iterable[str]:
    for spec in pkgutil.iter_modules(
            resolve_name(package_main_name).__spec__.submodule_search_locations):
        if not spec.name.startswith("_"):
            yield spec.name


def _get_doc(
        package_main_name: str,
        subcommand_name: str
) -> Optional[str]:
    """
    Return documentation of that module
    """
    importlib.import_module(f'{package_main_name}.{subcommand_name}')
    i = resolve_name(f'{package_main_name}.{subcommand_name}')
    if hasattr(i, '__doc__'):
        return i.__doc__
    else:
        return None


def get_main_func_from_subcommand(
        package_main_name: str,
        subcommand_name: str
) -> Optional[Callable[[List[str]], int]]:
    """
    Return a subcommands' "main" function.
    """
    importlib.import_module(f'{package_main_name}.{subcommand_name}')
    i = resolve_name(f'{package_main_name}.{subcommand_name}')
    if hasattr(i, 'main') and inspect.isfunction(getattr(i, 'main')):
        return i.main
    else:
        return None


def get_argparser_from_subcommand(
        package_main_name: str,
        subcommand_name: str
) -> Optional[argparse.ArgumentParser]:
    """
    Return result of a subcommands' "create_parser" function.
    """
    importlib.import_module(f'{package_main_name}.{subcommand_name}')
    i = resolve_name(f'{package_main_name}.{subcommand_name}')
    if hasattr(i, 'create_parser') and inspect.isfunction(getattr(i, 'create_parser')):
        return i.create_parser()
    else:
        return None

def lscmd(
        package_main_name: str,
        valid_subcommand_names: Iterable[str]
):
    _lh.info("Listing modules...")
    for item in valid_subcommand_names:
        doc = _get_doc(package_main_name, item)
        if doc is None:
            doc = NONE_DOC
        else:
            doc_splitlines = doc.splitlines()
            if not doc_splitlines:
                doc = NONE_DOC
            else:
                while len(doc_splitlines) > 0:
                    potential_doc = doc_splitlines[0].strip()
                    if potential_doc == "":
                        doc_splitlines.pop(0)
                    else:
                        doc = potential_doc
                        break
                else:
                    doc = NONE_DOC
        if doc.find("--") != -1:
            doc = doc.split("--")[1].strip()

        print(f"{item} -- {doc}")
    sys.exit(0)


class _ParsedArgs:
    input_subcommand_name: str = ""
    have_help: bool = False
    have_version: bool = False
    parsed_args: List[str] = []


def _parse_args(
        args: List[str]
) -> _ParsedArgs:
    parsed_args = _ParsedArgs()
    i = 0
    while i < len(args):
        name = args[i]
        if name in ('--help', '-h'):
            parsed_args.have_help = True
        elif name in ('--version', '-v'):
            parsed_args.have_version = True
        elif not name.startswith('-') and parsed_args.input_subcommand_name == "":
            parsed_args.input_subcommand_name = name
            args.pop(i)
        i += 1
    parsed_args.parsed_args = args
    return parsed_args


def _format_help_info(package_main_name: str) -> str:
    return f"""
This is frontend of `{package_main_name.split('.')[0].strip()}` provided by `commonutils.libfrontend`.

SYNOPSYS: {sys.argv[0]} [[SUBCOMMAND] [ARGS_OF SUBCOMMAND] ...] [-h|--help] [-v|--version]

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


def setup_frontend(
        package_main_name: str,
        one_line_description: str,
        version: str,
        help_info: Optional[str] = None,
        subcommand_help: str = "Use 'lscmd' to list all valid subcommands.",
        use_root_logger: bool = True,
        default_log_filename: str = "log.log"
):
    _lh.info(f'{one_line_description} ver. {version}')
    _lh.info(f'Called by: {" ".join(sys.argv)}')
    parsed_args = _parse_args(sys.argv[1:])
    if use_root_logger:
        log_filename = os.environ.get("LOG_FILE_NAME", default_log_filename)
        file_handler = logging.FileHandler(filename=log_filename)
        file_handler.setLevel(logger_helper.TRACE)
        file_handler.setFormatter(logger_helper.get_formatter(logger_helper.TRACE))
        logging.root.addHandler(file_handler)

    valid_subcommand_names = get_subcommands(package_main_name)
    if parsed_args.input_subcommand_name == "lscmd":
        lscmd(package_main_name, valid_subcommand_names)
    elif parsed_args.input_subcommand_name == "":
        if parsed_args.have_help:
            if help_info is None:
                help_info = _format_help_info(package_main_name)
            print(help_info)
            sys.exit(0)
        elif parsed_args.have_version:
            print(version)
            sys.exit(0)
        else:
            _lh.exception(f"Subcommand name not set! {subcommand_help}")
            sys.exit(1)
    elif parsed_args.input_subcommand_name in valid_subcommand_names:
        main_fnc = get_main_func_from_subcommand(package_main_name=package_main_name,
                                                 subcommand_name=parsed_args.input_subcommand_name)
        if main_fnc is not None:
            sys.exit(main_fnc(parsed_args.parsed_args))
        else:
            _lh.exception(f"Subcommand '{parsed_args.input_subcommand_name}' not found! {subcommand_help}")
            sys.exit(1)
    else:
        _lh.exception(f"Subcommand '{parsed_args.input_subcommand_name}' not found! {subcommand_help}")
        sys.exit(1)
