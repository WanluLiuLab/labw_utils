"""
labw_utils.devutils.sphinx_helper -- Helpers of Sphinx documentation system
"""
__all__ = (
    "convert_ipynb_to_myst",
)

import glob
import os
import re
from collections import defaultdict
from typing import List, Callable, Optional

import jupytext
import nbformat as nbf
import tomli

from labw_utils.commonutils import libfrontend
from labw_utils.commonutils.io.file_system import should_regenerate
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)

SHOUT_LINK_REGEX = re.compile(r"^# SHOUT LINK: (.+)$")

SHOUT_SEARCH_DICT = {
    "# RMCELL": "remove-cell",  # Remove the whole cell
    "# SKIP": "skip-execution",  # Remove the whole cell
    "# RMIN": "remove-input",  # Remove only the input
    "# RMOUT": "remove-output",  # Remove only the output
    "# HIDEIN": "hide-input"  # Hide the input w/ a button to show
}


def shell_filter(nb: nbf.NotebookNode) -> nbf.NotebookNode:
    for cell in nb.cells:
        if cell["cell_type"] != "code":
            continue
        cell_tags = cell.get('metadata', {}).get('tags', [])
        for key, val in SHOUT_SEARCH_DICT.items():
            if key in cell['source']:
                if val not in cell_tags:
                    cell_tags.append(val)
        if len(cell_tags) > 0:
            cell['metadata']['tags'] = cell_tags
    return nb


def convert_ipynb_to_myst(
        source_dir: str,
        hooks: Optional[List[Callable[[nbf.NotebookNode], nbf.NotebookNode]]] = None
):
    if hooks is None:
        hooks = []
    for fn in glob.glob(os.path.join(source_dir, "**", "*.ipynb"), recursive=True):
        fn = os.path.abspath(fn)
        if "ipynb_checkpoints" in fn or "_build" in fn:
            continue
        dst_fn = fn + ".md"
        if should_regenerate(fn, dst_fn):
            _lh.info(f"CONVERT {fn} -> {dst_fn}: START")
            try:
                nb = jupytext.read(fn, fmt="notebook")
                for hook in hooks:
                    nb = hook(nb)
                jupytext.write(
                    nb=nb,
                    fp=dst_fn,
                    fmt="md:myst"
                )
            except Exception as e:
                _lh.warning(f"CONVERT {fn} -> {dst_fn}: ERR", exc_info=e)
            _lh.info(f"CONVERT {fn} -> {dst_fn}: FIN")
        else:
            _lh.info(f"CONVERT {fn} -> {dst_fn}: REFUSE TO OVERWRITE NEWER FILE")


def generate_cli_docs(
        config_toml_file_path: str,
        dest_dir_path: str
):
    os.makedirs(dest_dir_path, exist_ok=True)
    with open(config_toml_file_path, "rb") as toml_reader:
        config_toml = tomli.load(toml_reader)

    arg_parsers = defaultdict(lambda: [])
    for main_module in config_toml["names"]:
        for subcommand in libfrontend.get_subcommands(main_module):
            parser = libfrontend.get_argparser_from_subcommand(main_module, subcommand)
            this_help_path = os.path.join(dest_dir_path, f"{main_module}.{subcommand}.txt")
            if parser is not None:
                with open(this_help_path, "w") as writer:
                    writer.write(parser.format_help())
                arg_parsers[main_module].append(subcommand)
            else:
                doc = libfrontend.get_doc_from_subcommand(main_module, subcommand)
                if doc is None:
                    continue
                else:
                    # doc_sio = io.StringIO(doc)
                    with open(this_help_path, "w") as writer:
                        writer.write(doc)
                    arg_parsers[main_module].append(subcommand)

    with open(os.path.join(dest_dir_path, "index.md"), "w") as index_writer:
        index_writer.write("# Command-Line Interfaces\n\n")
        for main_module, subcommands in arg_parsers.items():
            main_module_correct_name = main_module.replace("._main", "").replace(".main", "")
            index_writer.write(f"## `{main_module_correct_name}`\n\n")
            for subcommand in subcommands:
                index_writer.write(f"### `{subcommand}`\n\n")
                index_writer.write(
                    "```{literalinclude} " + f"{main_module}.{subcommand}.txt" + "\n:language: text\n```\n\n"
                )