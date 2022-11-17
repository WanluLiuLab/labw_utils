# On Standardization (Sorting) of GTF

## Introduction

[HUGO Gene Nomenclature Committee (HGNC)](https://www.genenames.org/)

Entrez

Ensembl format

Uniprot format

## External sorting algorithms

`bedtools`: Sorted by chromosome (alphabetical or self-provided) and coordinate (regardless of strand). Exon number unchanged.

`gffread`: Sorted by chromosome (alphabetical or self-provided) and coordinate (regardless of strand). Exon number depleted.

```{toctree}
:caption: 'Contents:'
:glob:
:maxdepth: 2

_gtf_sort/*
```
