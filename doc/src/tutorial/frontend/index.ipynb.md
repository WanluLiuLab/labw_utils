---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.5
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Frontend Tutorial

This tutorial provides a comprehensive introduction on frontends of `labw_utils`.

This tutorial assumes an `all` installation of optional extras or a full development installation. See [](../installation/index) for more details.

**How to read this documentation**: Code block without leading `!` are Python code blocks. For example:

```{code-cell}
import os; os.name
```

Code block with leading `!` are Shell code blocks. For example:

```{code-cell}
!ls -lFh | grep ipynb
```

## Preparation

```{code-cell}
import labw_utils

print(f"labw_utils: {labw_utils.__version__}")
```

Download data. Following would download an _C. Elegans_ TGS (ONT GridION) RNA-Seq data from ENA accession [ERR3245471](https://www.ebi.ac.uk/ena/browser/view/ERR3245471) (article {cite}`Roach2020`).

```{code-cell}
!if [ !-f L4_rep2.fastq.gz ]; then \
    wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3245471/L4_rep2.fastq.gz; \
else \
    echo "L4_rep2.fastq.gz already exists!" \
fi
```

## `labw_utils.bioutils`

The frontends underneath this package provides bioinformatics utilities. Use installation with `bioutils` extra to get best experience.

+++

### `labw_utils.bioutils lscmd`

This frontend does nothing but lists all available sub-commands and provides one-line description.

```{warning}
It would try to import every frontend so if you perform a default installation, frontends requiring optional extras like `get_exonic_depth` will NOT be shown.
```

```{code-cell}
!python -m labw_utils.bioutils lscmd
```

### `describe_fasta_by_binning`

Work in progress -- do not use.

+++

### `describe_fastq`

This is a light-weighted [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) implemented in pure Python. It supports NGS and TGS reads but without detection of repetitive sequences and adapters.

```{code-cell}
!
```
