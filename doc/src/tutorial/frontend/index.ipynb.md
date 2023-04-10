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
---
vscode:
  languageId: python
---
import os; os.name
```

Code block with leading `!` are Shell code blocks. For example:

```{code-cell}
---
vscode:
  languageId: python
---
!ls -lFh | grep ipynb
```

## Preparation

```{code-cell}
---
vscode:
  languageId: python
---
import labw_utils

print(f"labw_utils: {labw_utils.__version__}")
```

Download data. Following would download an _C. Elegans_ TGS (ONT GridION, R9.4 Pore) direct RNA-Seq data from ENA accession [ERR3245471](https://www.ebi.ac.uk/ena/browser/view/ERR3245471) (article {cite}`Roach2020`) and align it to UCSC `ce11` reference genome.

```{code-cell}
---
tags: [skip-execution]
vscode:
  languageId: python
---
# SKIP
!wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3245471/L4_rep2.fastq.gz
!gunzip L4_rep2.fastq.gz
!minimap2 -a -x splice ce11.fa L4_rep2.fastq | samtools sort -o L4_rep2.bam
!samtools index L4_rep2.bam
```

```{code-cell}
---
tags: [remove-input]
vscode:
  languageId: python
---
# RMIN
%cat preparation.log
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
---
vscode:
  languageId: python
---
!python -m labw_utils.bioutils lscmd
```

### `describe_fasta_by_binning`

Work in progress -- do not use.

+++

### `describe_fastq`

This is a light-weighted [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) implemented in pure Python. It supports NGS and TGS reads but without detection of repetitive sequences and adapters.

Help message:

```{code-cell}
---
vscode:
  languageId: python
---
!python -m labw_utils.bioutils describe_fastq
```

Following is an example of describing `L4_rep2.fastq`, the file we just downloaded.

```{code-cell}
---
tags: [skip-execution]
vscode:
  languageId: python
---
# SKIP
!python -m labw_utils.bioutils describe_fastq L4_rep2.fastq
```

```{code-cell}
---
tags: [remove-input]
vscode:
  languageId: python
---
# RMIN
%cat describe_fastq.log
```

This generates:

- `L4_rep2.fastq.stats.d`, the directory where quality control files are located. They are:
  - `all.tsv`, the per-read quality control file, with following columns:
    - `SEQID`, the FASTQ read ID.
    - `GC`, per-read GC content in absolute value.
    - `LEN`, actual read length.
    - `MEANQUAL`, mean sequencing quality using Phread33 score.
  - `extension_stat.tsv`, mean per-base quality of bases on each read from 5' to 3', mainly for NGS.
    - `POS`, position of each base on every transcript from 5' to 3'. If not present, would be omitted.
    - `QUAL`, mean sequencing quality using Phread33 score.

Example using `all.tsv`:

```{code-cell}
---
vscode:
  languageId: python
---

```
