# `labw_utils` -- Utility Python functions & classes used in LabW

**Markdown compatibility guide** This file is written in [Myst-flavored Markdown](https://myst-parser.readthedocs.io/), and may show errors on the default landing page of PYPI or Git Hostings. You can correctly preview it on generated Sphinx documentation or [Visual Studio Code](https://code.visualstudio.com) with [ExecutableBookProject.myst-highlight](https://marketplace.visualstudio.com/items?itemName=ExecutableBookProject.myst-highlight) plugin.

---

`labw_utils` contains a series of Python functions and classes for biological and general purpose programming in LabW.

The code base is designed with following principles:

1. Pure-Python implemented.
2. Minimal dependency. Contents inside this module try to not depend on third-party packages like Pandas or Numpy.

## Installation

### Using pre-built Library from PYPI

You need Python interpreter (CPython implementation) >= 3.8 (recommended 3.8) and latest [`pip`](https://pip.pypa.io/) to install this software from [PYPI](https://pypi.org). Command:

```shell
pip install labw_utils
```

You are recommended to use this application inside a virtual environment like [`venv`](https://docs.python.org/3/library/venv.html), [`virtualenv`](https://virtualenv.pypa.io), [`pipenv`](https://pipenv.pypa.io), [`conda`](https://conda.io) or [`poetry`](https://python-poetry.org).

### Build from Source

You need Python interpreter (CPython implementation) >= 3.8, latest PYPA [`build`](https://pypa-build.readthedocs.io) and latest [`setuptools`](https://setuptools.pypa.io/) to build this software. You are recommended to build the software in a virtual environment provided by [`virtualenv`](https://virtualenv.pypa.io), etc.

Build the software using:

```shell
rm -fr dist # Remove build of previous versions
python3 -m build
pip install dist/labw_utils-*.whl
```
