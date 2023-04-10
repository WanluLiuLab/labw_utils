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

This tutorial assumes basic Python and Shell scripting knowledge and understanding of underlying biological concepts like Gene, Transcript, Reference Genome, etc.

+++

**How to read this documentation**: Code block without leading `!` are Python code blocks. For example:

```{code-cell} ipython3
import time; print(time.asctime())
```

Code block with leading `!` are Shell code blocks. For example:

```{code-cell} ipython3
!ls -lFh | grep ipynb
```

## Preparation

```{code-cell} ipython3
# Development Block which can be safely ignored.
# %load_ext autoreload
# %autoreload 2
# %matplotlib inline

# pylint: disable=wrong-import-position, line-too-long, missing-module-docstring

import os
import sys

try:
    THIS_DIR_PATH = os.path.abspath(globals()["_dh"][0])
except KeyError:
    THIS_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
NEW_PYTHON_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(THIS_DIR_PATH)))), "src")
sys.path.insert(0, NEW_PYTHON_PATH)
os.environ["PYTHONPATH"] = os.pathsep.join((NEW_PYTHON_PATH, os.environ.get("PYTHONPATH", "")))
```

```{code-cell} ipython3
import labw_utils
import pandas as pd # Pandas for reading/writing relational data
import seaborn as sns # seaborn for simple plotting
import matplotlib.pyplot as plt # Another plotting library
import matplotlib

print(f"labw_utils: {labw_utils.__version__}")
print(f"pandas: {pd.__version__}")
print(f"seaborn: {sns.__version__}")
print(f"matplotlib: {matplotlib.__version__}")
```

Download data. Following would download an _C. Elegans_ TGS (ONT GridION, R9.4 Pore) direct RNA-Seq data from ENA accession [ERR3245471](https://www.ebi.ac.uk/ena/browser/view/ERR3245471) (article {cite}`Roach2020`) and align it to UCSC `ce11` reference genome.

```{code-cell} ipython3
:tags: [skip-execution]

# SKIP
!wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3245471/L4_rep2.fastq.gz
!gunzip L4_rep2.fastq.gz
!minimap2 -a -x splice ce11.fa L4_rep2.fastq | samtools sort -o L4_rep2.bam
!samtools index L4_rep2.bam
```

```{code-cell} ipython3
:tags: [remove-input]

# RMIN
%cat preparation.log
```

## `labw_utils.bioutils`

The frontends underneath this package provides bioinformatics utilities. Use installation with `bioutils` extra to get best experience.

+++

### `lscmd`

This frontend does nothing but lists all available sub-commands and provides one-line description.

```{warning}
It would try to import every frontend so if you perform a default installation, frontends requiring optional extras like `get_exonic_depth` will NOT be shown.
```

```{code-cell} ipython3
!python -m labw_utils.bioutils lscmd
```

### `describe_fasta_by_binning`

Work in progress -- do not use.

+++

### `describe_fastq`

This is a light-weighted [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) implemented in pure Python. It supports NGS and TGS reads but without detection of repetitive sequences and adapters.

Help message:

```{code-cell} ipython3
!python -m labw_utils.bioutils describe_fastq
```

Following is an example of describing `L4_rep2.fastq`, the file we just downloaded.

```{code-cell} ipython3
:tags: [skip-execution]

# SKIP
!python -m labw_utils.bioutils describe_fastq L4_rep2.fastq
```

```{code-cell} ipython3
:tags: [remove-input]

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

Following is an example of using `all.tsv` to determine base sequencing metrics.

```{code-cell} ipython3
fq_stats_all = pd.read_table("L4_rep2.fastq.stats.d/all.tsv")
fq_stats_all.head()
```

```{code-cell} ipython3
fig, axs = plt.subplots(3, 1)

sns.histplot(fq_stats_all, x="GC", ax=axs[0])
sns.histplot(fq_stats_all, x="LEN", ax=axs[1])
sns.histplot(fq_stats_all, x="MEANQUAL", ax=axs[2])
```

### `describe_gtf_by_binning`

Work in progress -- do not use.

+++

### `describe_gtf`

This shows basic QC metrics on GTF files that can be parsed into a Gene-Isoform-Exon three-tier structure.

TODO

+++

### `filter_gtf_by_attribute`

This is like the [`grep`](https://www.gnu.org/software/grep) utility for GTF. You are required to input an attribute name and a list of possible values and this function would filter GTF by that attribute.

````{note}
Attributes are the trailing `; ` separated key-value pair in a GTF record. For example, the GTF record 

```text
chrI	ncbiRefSeq	transcript	3747	3909	.	-	.	gene_id "Y74C9A.6"; transcript_id "NR_001477.2";  gene_name "Y74C9A.6";
```

have 3 attributes `gene_id` `transcript_id` `gene_name` with value `Y74C9A.6` `NR_001477.2` `Y74C9A.6`.

Fields like `chrI` or `3747` is more commonly refered to as "Required fields".
````

With following help message:

```{code-cell} ipython3
!python -m labw_utils.bioutils filter_gtf_by_attribute --help
```

For example, we wish to get all genomic features in [CED](https://wormbase.org/resources/gene_class/ced) gene class. It can be filtered by matching attribute `gene_id` with regular expression `^ced-.*$`. We would write a file named `filter.regex` with following contents:

```text
^ced-.*$
```

and perform filtering using:

```{code-cell} ipython3
:tags: [skip-execution]

# SKIP
!python -m labw_utils.bioutils filter_gtf_by_attribute \
    -g ce11.ncbiRefSeq.gtf \
    --attribute_name gene_id \
    --attribute_values filter.regex \
    --out ce11.filtered.gtf \
    --regex
```

```{code-cell} ipython3
:tags: [remove-input]

# RMIN
%cat filter_gtf_by_attribute.log
```

which generates `ce11.filtered.gtf` with CED gene class only.

+++

### `generate_fake_fasta`

Work in progress -- do not use.

+++

### `get_exonic_depth`

This script calculates sequencing depth of RNA-Seq reads. It would calculate number of bases of all primiarily mapped reads ovr mappable (exonic) region of provided genomic annotation.

```{warning}
This is designed for RNA-Seq only -- For Whole-Genome Sequencing (WGS) or Whole-Exon Sequencing (WES/WXS), please divide the number of bases by genome \& exon length.
```

```{warning}
This module does NOT take care of alt contigs, decoy sequences, patch fixes, unlocalized/unplaced scaffolds, etc. They would be regarded as individual chromosomes.
```

Help message:

```{code-cell} ipython3
!python -m labw_utils.bioutils get_exonic_depth --help
```

```{code-cell} ipython3
:tags: [skip-execution]

# SKIP
!python -m labw_utils.bioutils get_exonic_depth \
    -s L4_rep2.bam \
    -g ce11.ncbiRefSeq.gtf
```

```{code-cell} ipython3
:tags: [remove-input]

# RMIN
%cat get_exonic_depth.log
```

This indicate that the input alignment file, `L4_rep2.bam`, have 150,833,898 primiarily mapped base over 26,383,716 exonic bases on reference GTF. It have a sequencing depth of 5.72.

+++

### `normalize_gtf`

This scripts normalizes GTF into standard form.

TODO

+++

### `sample_transcript`

This script randomly samples isoforms in a GTF by percentage.

The help message is as follows:

```{code-cell} ipython3
!python -m labw_utils.bioutils sample_transcript --help
```

### `split_fasta`

This module splite FASTA with multiple contig into FASTA files with one file per contig.

The help message is as follows:

```{code-cell} ipython3
!python -m labw_utils.bioutils split_fasta --help
```

### `transcribe`

This step would transcribe the input genome GTF and genome FASTA into **stranded** transcriptome FASTA. It is designed to be general-purposed, i.e., can be applied on any matching GTF and FASTA. It should generate similar output with `bedtools getfasta -nameOnly -s -fi [FASTA] -bed [GTF] > [OUT]`

```{note}
Although this software can be used to generate reference cDNAs for software like Salmon, there are differences between transcribed cDNA and Ensembl-provided cDNA. Ensembl-provided cDNA does not include small features like lncRNA, while YASIM transcribed cDNA includes all transcripts inside provided GTF.
```

The help message is as follows:

```{code-cell} ipython3
!python -m labw_utils.bioutils transcribe --help
```

Example:

```{code-cell} ipython3
:tags: [skip-execution]

# SKIP
!python -m labw_utils.bioutils transcribe \
    -f ce11.fa \
    -g ce11.filtered.gtf \
    --no_write_single_transcript \
    -o ce11_trans_filtered.fa
```

```{code-cell} ipython3
:tags: [remove-input]

# RMIN
%cat transcribe.log
```

Generates:

- `ce11_trans_filtered.fa`, the generated cDNA sequence FASTA.
- `ce11_trans_filtered.fa.stats`, a TSV file with following columns:
  - `TRANSCRIPT_ID`, the `transcript_id` field in GTF.
  - `GENE_ID`, the `gene_id` field in GTF.
  - `SEQNAME`, chromosome or contig name.
  - `START`, the `start` field in GTF, 1-based inclusive.
  - `END`, the `end` field in GTF, 1-based inclusive.
  - `STRAND`, the `strand` field in GTF.
  - `ABSOLUTE_LENGTH`, is `START` - `END` + 1.
  - `TRANSCRIBED_LENGTH`, length of the cDNA without introns and UTRs.
  - `GC`, GC content of the cDNA in percentage.

If `--no_write_single_transcript` option was turned off, following additional directory will be generated:

- `ce11_trans_filtered.fa.d`, the directory where every cDNA is stored as separate FASTA.

+++

## `ysjs`

YSJS -- YuZJ Simple Job System, a [SLRUM](https://slurm.schedmd.com/overview.html) that runs on single machine.

Work in progress -- do not use.

+++

## `ysjsd`

The server side of YSJS.

Work in progress -- do not use.
