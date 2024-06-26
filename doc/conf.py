"""
Configuration file for the Sphinx documentation builder.
"""

# pylint: disable=wrong-import-position, invalid-name

import os

from docutils.parsers.null import Parser as NullParser
from sphinx.application import Sphinx

import labw_utils
from labw_utils.stdlib.cpy311 import tomllib

os.environ["LABW_UTILS_SPHINX_BUILD"] = "1"  # Disable chronolog and others.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(THIS_DIR)


def setup(app: Sphinx):
    app.add_source_parser(NullParser)


# -- Project information -----------------------------------------------------

with open(os.path.join(ROOT_DIR, "pyproject.toml"), "rb") as reader:
    parsed_pyproject = tomllib.load(reader)

project = parsed_pyproject["project"]["name"]
author = "&".join([author["name"] for author in parsed_pyproject["project"]["authors"]])
copyright_string = f"2022-2023, {author}"
release = labw_utils.__version__

# -- General configuration ---------------------------------------------------

html_theme = parsed_pyproject["tool"]["sphinx"]["html_theme"]
extensions = parsed_pyproject["tool"]["sphinx"]["extensions"]

myst_enable_extensions = ["deflist", "colon_fence"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", ".virtualenv/**"]

# html_static_path = ['_static']

# Source code suffixes
source_suffix = {".rst": "restructuredtext", ".md": "myst-nb", ".ipynb": "null"}
nb_custom_formats = {
    ".ipynb.py": ["jupytext.reads", {"fmt": "py:percent"}],
    ".ipynb.md": ["jupytext.reads", {"fmt": "md:myst"}],
}

# Autodoc settings
autodoc_default_options = {
    "special-members": "__init__",
}
autodoc_class_signature = "separated"
autodoc_member_order = "bysource"
autodoc_typehints = "description"

# Intersphinx settings
intersphinx_pkgs = [
    "python",
    "joblib",
    "torch",
    "psutil",
    "pandas",
    "numpy",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3.8", ("intersphinx/python.inv", None)),
    "joblib": ("https://joblib.readthedocs.io/en/latest", ("intersphinx/joblib.inv", None)),
    "torch": ("https://pytorch.org/docs/stable", ("intersphinx/torch.inv", None)),
    "psutil": ("https://psutil.readthedocs.io/en/latest", ("intersphinx/psutil.inv", None)),
    "pandas": ("https://pandas.pydata.org/docs/", ("intersphinx/pandas.inv", None)),
    "numpy": ("https://numpy.org/doc/stable/", ("intersphinx/numpy.inv", None)),
}
manpages_url = "https://manpages.debian.org/{path}"

# myst-nb settings
nb_execution_timeout = 1200
if os.getenv("NO_NB") is not None:
    nb_execution_mode = "off"
else:
    nb_execution_mode = "cache"
nb_merge_streams = True

# BibTeX setting
bibtex_bibfiles = ["refs.bibtex.bib"]
bibtex_default_style = "plain"
bibtex_reference_style = "super"

# sphinx-hoverxref from <https://github.com/scverse/scvi-tools/blob/main/docs/conf.py>
hoverx_default_type = "tooltip"
hoverxref_domains = ["ref", "class", "func", "meth", "attr", "exc", "data", "mod", "py", "obj"]
hoverxref_role_types = dict.fromkeys(
    ["ref", "class", "func", "meth", "attr", "exc", "data", "mod", "py", "obj"],
    "tooltip",
)
hoverxref_intersphinx = intersphinx_pkgs
