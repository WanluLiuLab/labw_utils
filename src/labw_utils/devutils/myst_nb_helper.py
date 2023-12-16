"""
``labw_utils.devutils.myst_nb_helper`` -- Helpers of Sphinx documentation system

TODO: This module is largely unfinished.

.. warning::
    UNSTABLE -- SUBJECT TO CHANGE!
"""

__all__ = ("convert_ipynb_to_myst", "shell_filter", "generate_cli_docs")

import glob
import os
import re
from collections import defaultdict

from labw_utils import UnmetDependenciesError, PackageSpecs, PackageSpec
from labw_utils.commonutils import libfrontend
from labw_utils.commonutils.lwio.file_system import should_regenerate
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.stdlib.cpy311 import tomllib
from labw_utils.typing_importer import Literal, Optional, Callable, List

PackageSpecs.add(
    PackageSpec(
        name="jupytext",
        conda_name="jupytext",
        conda_channel="conda-forge",
        pypi_name="jupytext",
    )
)
PackageSpecs.add(
    PackageSpec(
        name="nbformat",
        conda_name="nbformat",
        conda_channel="conda-forge",
        pypi_name="nbformat",
    )
)

try:
    import pytest

    jupytext = pytest.importorskip("jupytext")
    nbf = pytest.importorskip("nbformat")
except ImportError:
    pytest = None
    try:
        import jupytext
    except ImportError as e:
        raise UnmetDependenciesError("jupytext") from e
    try:
        import nbformat as nbf
    except ImportError as e:
        raise UnmetDependenciesError("nbformat") from e

_lh = get_logger(__name__)

SHOUT_LINK_REGEX = re.compile(r"^# SHOUT LINK: (.+)$")

SHOUT_SEARCH_DICT = {
    "# RMCELL\n": "remove-cell",  # Remove the whole cell
    "# SKIP\n": "skip-execution",  # Remove the whole cell
    "# RMIN\n": "remove-input",  # Remove only the input
    "# RMOUT\n": "remove-output",  # Remove only the output
    "# HIDEIN\n": "hide-input",  # Hide the input w/ a button to show
}


def shell_filter(nb: nbf.NotebookNode) -> nbf.NotebookNode:
    for cell in nb.cells:
        if cell["cell_type"] != "code":
            continue
        cell_tags = cell.get("metadata", {}).get("tags", [])
        for key, val in SHOUT_SEARCH_DICT.items():
            if key in cell["source"]:
                if val not in cell_tags:
                    cell_tags.append(val)
                cell["source"] = cell["source"].replace(key, "")
        if len(cell_tags) > 0:
            cell["metadata"]["tags"] = cell_tags
    return nb


def convert_ipynb_to_myst(
    source_dir: str,
    hooks: Optional[List[Callable[[nbf.NotebookNode], nbf.NotebookNode]]] = None,
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
                jupytext.write(nb=nb, fp=dst_fn, fmt="md:myst")
            except Exception as _e:
                _lh.warning(f"CONVERT {fn} -> {dst_fn}: ERR", exc_info=_e)
            _lh.info(f"CONVERT {fn} -> {dst_fn}: FIN")
        else:
            _lh.info(f"CONVERT {fn} -> {dst_fn}: REFUSE TO OVERWRITE NEWER FILE")


def generate_cli_docs(
    config_toml_file_path: str,
    dest_dir_path: str,
    format: Literal["txt", "myst.md"] = "txt",
):
    os.makedirs(dest_dir_path, exist_ok=True)
    with open(config_toml_file_path, "rb") as toml_reader:
        config_toml = tomllib.load(toml_reader)

    arg_parsers = defaultdict(lambda: [])
    for main_module in config_toml["names"]:
        for subcommand in libfrontend.get_subcommands(main_module):
            print(f"Generating CLI docs for {subcommand}")
            parser = libfrontend.get_argparser_from_subcommand(main_module, subcommand)
            this_help_path = os.path.join(dest_dir_path, f"{main_module}.{subcommand}.{format}")
            if parser is not None:
                with open(this_help_path, "w") as writer:
                    if format == "myst.md":
                        raise NotImplementedError
                        # writer.write(parser.to_markdown())
                    elif format == "txt":
                        writer.write(parser.format_help())
                arg_parsers[main_module].append(subcommand)
            else:
                doc = libfrontend.get_doc_from_subcommand(main_module, subcommand)
                if doc is None:
                    continue
                else:
                    with open(this_help_path, "w") as writer:
                        writer.write(doc)
                    arg_parsers[main_module].append(subcommand)

    with open(os.path.join(dest_dir_path, "index.md"), "w") as index_writer:
        index_writer.write("# Command-Line Interfaces\n\n")
        for main_module, subcommands in arg_parsers.items():
            main_module_correct_name = main_module.replace("._main", "").replace(".main", "")
            index_writer.write(f"## `{main_module_correct_name}`\n\n")
            for subcommand in subcommands:
                index_writer.write(f"### `{main_module_correct_name}` `{subcommand}`\n\n")
                if format == "myst.md":
                    index_writer.write("```{include} " + f"{main_module}.{subcommand}.{format}" + "\n```\n\n")
                elif format == "txt":
                    index_writer.write(
                        "```{literalinclude} " + f"{main_module}.{subcommand}.{format}" + "\n:language: text\n```\n\n"
                    )
