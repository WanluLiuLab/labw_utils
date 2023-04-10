# Installing `labw_utils`

`labw_utils` is a package that can both be used as an executable and as a supporting library. Here provides a detailed guide on how it should be installed. In following text, `[VERSION]` should be replaced with your desired version, currently 1.0.1.

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

```shell
python -m labw_utils.bioutils lscmd
```

which should give an output like:

```text
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
    - pytorch
    - nvidia

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

You're then free to go! Your editors or IDEs should recognize `labw_utils` and give you type hints \& hover docs afterwards.

## Dependencies

The default installation requires following packages:

- `tqdm`, for a user-friendly progress bar.
  - Although not recommended, this package can be uninstalled. If so, `labw_utils` would use a bundled shabby process bar implementation.
- `numpy`, for performing numerical operations. This package is used in most places so uninstallation is not recommended.
- `tomli` and `tomli_w`, for parsing and serialization from \& to TOML formats.
  - These packages can be uninstalled if you do not use {py:mod}`labw_utils.commonutils.serializer.toml` module.
- `pandas`, for reading and writing relational data formats. This package is used in most places so uninstallation is not recommended.
- `tomli` and `tomli_w`, for parsing and serialization from \& to TOML formats.
- `joblib`, for embarrasing parallelization of small tasks.

### Non-Default Installation with Optional Extras

Although the default packing stratergy can satisfy most biological uses, you may see errors like:

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

in place where default stratergy fails. You can use pip's "optional extra" function to install other dependencies. For example, to perform installation with `bioutils` extra, you should:

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

- `pysam`, for parsing Sequence Alignment Format (SAM)/Binary Alignment Format (BAM) files.
  - Notice: This package is not installable under Microsoft Windows.
  - Used in `get_exonic_depth` executable.

#### `mlutils`

The `mlutils` installation come with following additional package:

- `torch`, for pyTorch integration. With this module:
  - The {py:mod}`labw_utils.mlutils.io_helper` would work with pyTorch Tensor serialization format.
  - The {py:mod}`labw_utils.mlutils.ndarray_helper` would work with {py:class}`torch.Tensor`.
  - The {py:mod}`labw_utils.mlutils.torch_layers` would function.

#### `ysjs`

The `ysjs` installation would allow user to use YSJS Client.

#### `ysjsd`

The `ysjsd` installation would allow user to use YSJSD Server.

#### `appenders`

The `appenders` installation come with following additional package:

- `pyarrow`, to use Apache Parquet as output format in `describe_fasta_by_binning` and `describe_gtf_by_binning` executables.
- `fastparquet`, to use Apache Parquet format in {py:mod}`labw_utils.commonutils.appender`.
- `pytables`, to use HDF5 format in {py:mod}`labw_utils.commonutils.appender`.

#### `all`

The `all` installation installs all above.