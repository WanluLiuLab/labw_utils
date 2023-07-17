# Installing `labw_utils`

`labw_utils` is a package that can both be used as an executable and as a supporting library. Here provides a detailed guide on how it should be installed. In following text, `[VERSION]` should be replaced with your desired version, currently 1.0.1.

```{warning}
Copy-and-paste code on a website to your terminal is dangerous. See following StackOverflow pages for more details on prevention:

- [web browser - Simple way to safely paste text from website into terminal - Information Security Stack Exchange](https://security.stackexchange.com/questions/249586/simple-way-to-safely-paste-text-from-website-into-terminal)
- [malware - What is the risk of copy and pasting Linux commands from a website? How can some commands be invisible? - Information Security Stack Exchange](https://security.stackexchange.com/questions/249586/simple-way-to-safely-paste-text-from-website-into-terminal)
- [exploit - How can I protect myself from this kind of clipboard abuse? - Information Security Stack Exchange](https://security.stackexchange.com/questions/39118/how-can-i-protect-myself-from-this-kind-of-clipboard-abuse)
```

## Prerequisites

The `labw_utils` is implemented in Python, so a working Python intepreter is required. Following are several popular ways to get it installed:

::::{tab-set}

:::{tab-item} Conda
:sync: ts_k_conda

We assume you're on a machine with [Conda](https://conda.io) (or [Mamba](https://mamba.readthedocs.io/)) installed. If not, get it from [Here](https://docs.conda.io/en/latest/miniconda.html). On successful installation, you should be able to execute following command in your terminal:

```console
$ conda --version
conda 23.3.1
```

:::

:::{tab-item} pip
:sync: ts_k_pip

To use PYPA [`pip`](https://pip.pypa.io/), you should install Python interpreter first.

Firstly, make sure that you're using correct interpreter using:

```console
$ which python
/home/yuzj/conda/envs/labw_utils/bin/python
```

Check your Python version using:

```console
$ python --version
Python 3.8.15
```

If the command is failed or if the version is too low, get one from [Official Site](https://www.python.org) or your package management systems.

Then you may check whether `pip` is installed using:

```console
$ python -m pip --version
pip 23.0.1 from /home/yuzj/conda/envs/labw_utils/lib/python3.8/site-packages/pip (python 3.8)
```

If failed, you may either install `pip` using package manager (If you installed Python that way) or use script from <https://bootstrap.pypa.io/get-pip.py>.

```{warning}
Recommended to use `python -m pip` instead of `pip` directly since the first `pip` your Shell found may belongs to another Python intepreter.
```

:::

::::

## Installation as an Executable

This section is for those who uses executables inside {py:mod}`labw_utils.bioutils`, like `transcribe` otr `describe_fastq`.

::::{tab-set}

:::{tab-item} Conda
:sync: ts_k_conda

You are recommended to use `labw_utils` in a separate environment since it may pollute your `base` Conda environment. Use:

```shell
conda create -n labw_utils python=3.8
codna activate labw_utils
pip install labw_utils~=[VERSION]
```

:::

:::{tab-item} pip
:sync: ts_k_pip

You are recommended to use `labw_utils` in a separate environment since it may pollute your global or user `site-packages` (or `dist-packages`, as-is in Debian-based GNU/Linux). You may create virtual environment using [`venv`](https://docs.python.org/3/library/venv.html), [`virtualenv`](https://virtualenv.pypa.io), [`pipenv`](https://pipenv.pypa.io), [`conda`](https://conda.io) or [`poetry`](https://python-poetry.org). You are recommended to use `venv` as it is usually binding with your Python installation [^Debian]. After activating your virtual environment, use:

```shell
pip install labw_utils~=[VERSION]
```

[^Debian]: In Debian-based GNU/Linux, install [`python3-venv`](https://packages.debian.org/stable/python3-venv) before proceeding if you installed your Python using APT.
:::

::::

You may test whether your installation is successful using:

```console
$ python -m labw_utils.bioutils lscmd
2023-04-09 21:18:38,800	[INFO] labw_utils.bioutils -- Biological Utilities used in LabW projects ver. 1.0.1
2023-04-09 21:18:38,800	[INFO] Called by: /home/yuzj/Documents/labw_utils/src/labw_utils/bioutils/__main__.py lscmd
2023-04-09 21:18:38,800	[INFO] Listing modules...
2023-04-09 21:18:38,940	[INFO] Note: NumExpr detected 36 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
2023-04-09 21:18:38,940	[INFO] NumExpr defaulting to 8 threads.
2023-04-09 21:18:39,061	[WARNING] You are importing from gene_view_v0_1_x -- This GTF parser is NOT safe and will be deprecated.
2023-04-09 21:18:39,064	[WARNING] This module is not finished
describe_fasta_by_binning -- Describe statistics on FASTA by binning.
describe_fastq -- Lite Python-implemented fastqc.
describe_gtf -- Get statistics about GTF files that can be parsed into a Gene-Transcript-Exon Three-Tier Structure
describe_gtf_by_binning -- Describe number of features in GTF by binning.
generate_fake_fasta -- Create fake organism.
get_exonic_depth -- Get depth from RNA-Seq alignment files.
get_transcript -- Filter GTF records by a specific attributes
normalize_gtf -- Performs GTF normalization, etc.
sample_transcript -- Sample fraction of transcripts inside GTF file.
split_fasta -- Split input FASTA file into one-line FASTAs with one file per contig.
transcribe -- General-purposed stranded transcription, from reference genome to reference cDNA.
```

## Installation as a Library

This section is for those who uses functions or classes in `labw_utils` to develop their own software.

Suppose you're developing a software package named `foo` using Python 3.8 and wish to depend itself on `labw_utils`, you should do as follows:

::::{tab-set}

:::{tab-item} Conda
:sync: ts_k_conda

1. Create a file named `foo_dev.yml` as follows:

    ```yaml
    name: foo_dev
    channels:
    - nodefaults
    - conda-forge
    - anaconda
    - bioconda
    - main
    - free

    dependencies:
    # Core libraries
    - python=3.8
    - [OTHER CONDA DEPENDENCIES]

    - pip:
        - labw_utils~=[VERSION]
        - [OTHER PIP DEPENDENCIES]
    ```

2. Create this environment using `conda env create -f foo_dev.yml`.
3. In your editor or IDE, activate this environment.
:::

:::{tab-item} pip
:sync: ts_k_pip

1. Create a file named `requirements.txt` as follows:

    ```text
    labw_utils~=[VERSION]
    ```

2. Create virtual environment.
3. Activate the virtual environmtnt, and install package using `pip install -r requirements.txt`
4. In your editor or IDE, activate this environment.
:::

::::

You may test whether your installation is successful using:

```pycon
>>> import labw_utils
>>> labw_utils.__version__
1.0.1
```

## Dependencies

The default installation requires following packages:

- [`tqdm`](https://tqdm.github.io), for a user-friendly progress bar.
  - Although not recommended, this package can be uninstalled. If so, `labw_utils` would use a bundled shabby process bar implementation.
- [`numpy`](https://numpy.org), for performing numerical operations. This package is used in most places so uninstallation is not recommended.
- [`tomli-w`](https://pypi.org/project/tomli-w), for serialization to TOML formats.
  - These packages can be uninstalled if you do not use {py:mod}`labw_utils.commonutils.serializer.toml` module.
- [`pandas`](https://pandas.pydata.org), for reading and writing relational data formats. This package is used in most places so uninstallation is not recommended.
- [`joblib`](https://joblib.readthedocs.io), for embarrasing parallelization of small tasks.

### Non-Default Installation with Optional Extras

Although the default packing strategy can satisfy most biological uses, you may see errors like:

```pycon
>>> from labw_utils.mlutils import torch_layers
Traceback (most recent call last):
  File "/home/yuzj/Documents/labw_utils/src/labw_utils/mlutils/torch_layers.py", line 9, in <module>
    import torch
ModuleNotFoundError: No module named 'torch'

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/yuzj/Documents/labw_utils/src/labw_utils/mlutils/torch_layers.py", line 11, in <module>
    raise UnmetDependenciesError("torch") from e
labw_utils.UnmetDependenciesError: torch not installed; Use ``conda install -c pytorch pytorch``; Use ``pip install torch``
```

in place where default strategy fails. You can use pip's "optional extra" function to install other dependencies. For example, to perform installation with `bioutils` extra, you should:

```shell
pip install labw_utils[bioutils]~=[VERSION]
```

which also works in `requirements.txt`.

You can specify multiple extras. For example:

```shell
pip install labw_utils[bioutils,mlutils]~=[VERSION]
```

A list of extras are:

#### `bioutils`

The `bioutils` installation come with following additional packages:

- [`pysam`](https://pysam.readthedocs.io), for parsing Sequence Alignment Format (SAM)/Binary Alignment Format (BAM) files.
  - Used in `get_exonic_depth` executable.

```{warning}
This installation option is not installable under Microsoft Windows.
```

#### `mlutils`

The `mlutils` installation come with following additional package:

- [`torch`](https://pytorch.org), for pyTorch integration. With this module:
  - The {py:mod}`labw_utils.mlutils.io_helper` would work with pyTorch Tensor serialization format.
  - The {py:mod}`labw_utils.mlutils.ndarray_helper` would work with {py:class}`torch.Tensor`.
  - The {py:mod}`labw_utils.mlutils.torch_layers` would function.

#### `ysjs`

The `ysjs` installation would allow user to use YSJS Client. That insluces:

- [`requests`](https://docs.python-requests.org) for {py:mod}`libysjs.operation`

#### `ysjsd`

The `ysjsd` installation would allow user to use YSJSD Server. That insluces:

- [`flask`](https://flask.palletsprojects.com) for backend HTTP server.
- [`sqlalchemy`](https://www.sqlalchemy.org) for connecting to databases.
- [`psutil`](https://psutil.readthedocs.io/en/latest) for resource monitoring.
- [`gevent`](https://www.gevent.org) for CGI server.
- [`jinja2`](https://jinja.palletsprojects.com) for template rendering.

#### `appenders`

The `appenders` installation come with following additional package:

- [`pyarrow`](https://arrow.apache.org/docs/python), to use [Apache Parquet](https://parquet.apache.org) as output format in `describe_fasta_by_binning` and `describe_gtf_by_binning` executables.
- [`fastparquet`](https://fastparquet.readthedocs.io), to use Apache Parquet format in {py:mod}`labw_utils.commonutils.appender`.
- [`pytables`](https://www.pytables.org), to use HDF5 format in {py:mod}`labw_utils.commonutils.appender`.

#### `all`

The `all` installation installs all above.

## Alternate Python Implementation and Low Python Versions

Although using Python 3.8 is recommended, it is also possible to use most function of `labw_utils` on Python 3.7. You need to specify `--ignore-requires-python` flag in PIP.

Known limitations are:

1. The typing system would be completely a mess. The {py:class}`typing.Final` and {py:class}`typing.Literal` would not work as expected.
2. Version of dependent packages (i.e., Pandas, Numpy, etc.) would be low, making it impossible for you to enjoy advantages provided by more recent versions.

Python3 less than 3.6 and Python <= 3 is explicitly unsupported. Do not try on these Python implementations. [Python implementation](https://wiki.python.org/moin/PythonImplementations) except from [CPython](https://github.com/python/cpython) and [PyPy](https://www.pypy.org) are not recommended. Use at your own risk.

Following is a list of tested alternate Python implementations and versions using [Tox](https://tox.wiki/). Dependencies are set up using either `all` or default installation strategy using Conda.

| Implementation | Version | Python API | ALL                | DEFAULT            |
|----------------|---------|------------|--------------------|--------------------|
| CPython        | NA      | 3.6        | NOT TESTED         | NOT TESTED         |
| CPython        | 3.7.12  | 3.7        | PASS               | PASS               |
| CPython        | 3.8.16  | 3.8        | PASS               | PASS               |
| CPython        | 3.9.16  | 3.9        | PASS               | PASS               |
| CPython        | 3.10.10 | 3.10       | PASS               | PASS               |
| CPython        | 3.11.3  | 3.11       | UNMET DEPENDENCIES | PASS               |
| CPython        | NA      | 3.12       | NOT TESTED         | NOT TESTED         |
| PyPy           | NA      | 3.6        | NOT TESTED         | NOT TESTED         |
| PyPy           | 7.3.7   | 3.7        | PASS               | PASS               |
| PyPy           | 7.3.11  | 3.8        | UNMET DEPENDENCIES | FAILED             |
| PyPy           | NA      | 3.9        | UNMET DEPENDENCIES | UNMET DEPENDENCIES |
| GraalPy        |         | 3.8        | UNMET DEPENDENCIES | UNMET DEPENDENCIES |

Side-loaded acceleration techniques:

| Implementation | Version | Python API | ALL   | DEFAULT |
|----------------|---------|------------|-------|---------|
| Pyston         | 2.3.4   | 3.8        | PASS  | PASS    |
| Pyjion         | 1.2.6   | 3.10       | PASS  | FAILED  |

Python API supported by [Jython](https://www.jython.org) and [IronPython](https://ironpython.net) is too low, so not assessed.
