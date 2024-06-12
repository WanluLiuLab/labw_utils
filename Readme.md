# `labw_utils` -- Utility Python functions & classes used in LabW

**Markdown compatibility guide** This file is written in [Myst-flavored Markdown](https://myst-parser.readthedocs.io/), and may show errors on the default landing page of PYPI or Git Hosting. You can correctly preview it on generated Sphinx documentation or [Visual Studio Code](https://code.visualstudio.com) with [ExecutableBookProject.myst-highlight](https://marketplace.visualstudio.com/items?itemName=ExecutableBookProject.myst-highlight) plugin.

---

**Not finished -- Do not use.**

Badages:
[![Python version](https://img.shields.io/badge/python-3.8-blue)](https://www.python.org/downloads/)
[![PyPI - Version](https://img.shields.io/pypi/v/labw_utils)](https://pypi.org/project/labw_utils/)
[![GitHub Contributors](https://img.shields.io/github/contributors/WanluLiuLab/labw_utils)](https://github.com/WanluLiuLab/labw_utils)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![License](https://img.shields.io/badge/licence-MIT-blue)](https://mit-license.org/)

URLs: [PYPI](https://pypi.org/project/yasim/), [GitHub](https://github.com/WanluLiuLab/yasim).

`labw_utils` contains a series of [Python](http://www.python.org/) functions and classes for biological and general purpose programming in [LabW](https://labw.org/).

The code base is designed with the following principles:

1. Pure-Python implemented.
2. Minimal dependency.

## Installation

### Using pre-built Library from PYPI

You need Python interpreter (CPython implementation) >= 3.8 (recommended 3.9) and latest [`pip`](https://pip.pypa.io/) to install this software from [PYPI](https://pypi.org). Command:

```shell
pip install labw_utils[defaults]
```

You are recommended to use this application inside a virtual environment like [`venv`](https://docs.python.org/3/library/venv.html), [`virtualenv`](https://virtualenv.pypa.io), [`pipenv`](https://pipenv.pypa.io), [`conda`](https://conda.io) or [`poetry`](https://python-poetry.org).

### Build from Source

You need Python interpreter (CPython implementation) >= 3.8 (recommended 3.9), latest PYPA [`build`](https://pypa-build.readthedocs.io) and latest [`setuptools`](https://setuptools.pypa.io/) to build this software. You are recommended to build the software in a virtual environment provided by [`virtualenv`](https://virtualenv.pypa.io), etc.

Build the software using:

```shell
rm -fr dist # Remove build of previous versions
python3 -m build
pip install dist/labw_utils-*.whl
```
