# FASTA and Genome Assembly Behind FASTA

## Introduction

The choise of reference genome GTF and FASTA is important in analysis process. Here we would discuss differences among these FASTA, reason behind them and how we might choose the right version.

## Piror Knowledge on Genome Assembly

Before we have reference FASTA, we would assembly the genome into chromosome and other sequences. This is implemented using following way:

- Clone sequences, Whole Genome Sequencing (WGS) sequences or PCR fragments termed **components** are retrived from sequencer and submitted to NCBI GenBank {cite}`GRC_assembly_terminology`.
- Components are **joint** into **contigs**, which are contineous genomic regions {cite}`GRC_assembly_terminology`.
- Ordered gapped or ungapped contigs are assembled into **scaffolds**. At that time their location on genome is unknown {cite}`GRC_assembly_terminology`.
- Scaffolds are assembled into **chromosomes** {cite}`GRC_assembly_terminology`. There might be **unplaced** (aka., **random**. known to belong on a specific chromosome but with unknown order or orientation) or **unlocalized** (chromosome of origin unknown) scaffolds {cite}`CaetanoAnolles2022`.
- Chromosomes, together with unlocalized and unplaced contigs and scaffolds, are packed together and released as **assemblies** {cite}`GRC_assembly_terminology`.

There are different types of assemblies. They are:

Primiary Assembly
: The assembly from one man. This excludes presense of alternative loci. A loci may be sequenced one or zero time {cite}`GRC_assembly_terminology`.

Haploid Assembly
: Primiary assembly together with alternative loci. A loci may be sequenced zero, one or more time but no duplicated chromosome whould present {cite}`GRC_assembly_terminology`.

Diploid Assembly
: Diploid version of haploid assembly, should contain 2 set of chromosomes without alternative loci {cite}`GRC_assembly_terminology`.

Toplevel Assembly
: Assembly that contained only top-level sequences of overlapping regions. i.e. contigs are not included if they are part of scaffolds, scaffolds are not included if they are part of chromosomes {cite}`NCBI_assembly_help`.

Beside the chromosome, an assembly may also include following sequcnes:

Alternative loci
: aka., "partial chromosomes", "alternate alleles", and "alternate haplotypes", which provides an alternate representation of a locus found in a largely haploid assembly {cite}`GRC_assembly_terminology`.

FIX patch
: A patch that corrects errors in current sequencing result {cite}`GRC_assembly_terminology`.

NOVEL patch
: A patch that finds another new alternative loci {cite}`GRC_assembly_terminology`.

FIX and NOVEL patches are added in minor release (e.g., GRCh38.p13 -> GRCh38.p14) and would be eliminated in major release (e.g., GRCh37 -> GRCh38)  {cite}`GRC_assembly_terminology,NCBI_assembly_help`.

For human genomes, there are also some unique features, including:

Pseudo-Autosomal Regions (PAR) and Satelites
: There are also duplications on centromeric regions and DNA satelites on chromosomes 5, 14, 19, 21 \& 22 {cite}`CaetanoAnolles2022`.

Epstein-Barr Virus (EBV) and Decoy sequences
: Some people would include EBV and Decoy sequences to enhance mapping in SNP calling pipelines. EBV is a human endogenerous retrovirus (hERV) that infects B-cells and are commonly seen in sequencing results, while decoy sequences are a set of highly repetitive genome sequences missing in human genome references {cite}`Minikel2013,Li2017` that are considered as a sink for alignment of reads that are often present in sequencing samples {cite}`NCBI_analysis_set`.

Cambridge Reference Sequence (rCRS) for Mitochondrial
: rCRS {cite}`Anderson1981` is a popular mitochondrial sequence used for population genetics. The official GRCh37 comes with a mitochondrial sequence 2bp longer than rCRS while GRCh38 do not have this problem {cite}`Li2017`.

International Union of Biochemists (IUB) Codes
: Codes to represent non-AGCT nucleotides. They are (reproduced from <https://biocorp.ca/IUB.php>):

> | Code | Definition | Mnemonic       |
> | ---- | ---------- | -------------- |
> | A    | Adenine    | **A**          |
> | C    | Cytosine   | **C**          |
> | G    | Guanine    | **G**          |
> | T    | Thymine    | **T**          |
> | R    | AG         | pu**R**ine     |
> | Y    | CT         | p**Y**rimidine |
> | K    | GT         | **K**eto       |
> | M    | AC         | a**M**ino      |
> | S    | GC         | **S**trong     |
> | W    | AT         | **W**eak       |
> | B    | CGT        | Not A          |
> | D    | AGT        | Not C          |
> | H    | ACT        | Not G          |
> | V    | ACG        | Not T          |
> | N    | AGCT       | a**N**y        |

## Chromosome Name Conventions

NCBI RefSeq Accession
: The chromosome number is the accession used in NCBI RefSeq. They are:
    - Chromosome names are in NCBI RefSeq Complete genomic molecule names. That is: `NC_0000[ID].[VER]`. For example, [NC_000001.11](https://www.ncbi.nlm.nih.gov/nuccore/NC_000001.11) for chromosome 1 version 11. [`NC_000023`](https://www.ncbi.nlm.nih.gov/nuccore/NC_000023) [`NC_000024`](https://www.ncbi.nlm.nih.gov/nuccore/NC_000024) [`NC_012920`](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920) are for chromosome X, Y and mitochondrion {cite}`NCBI2013`.
    - Unplaced/unlocalized scaffolds of a chromosome are started with `NT_` (clone-based or WGS Genomic Contig or scaffold) after the chromosome.
    - Genomic patch and alternative loci are satrted with `NW_` (WGS Genomic Contig or scaffold) or `NT_` after all scaffolds {cite}`NCBI2013`.

GenBank Accession
: The GenBank Accession starts with `CM`. For example, [CM000663](https://www.ncbi.nlm.nih.gov/nuccore/CM000663) for chromosome 1.

Analysis Set Naming Convention
: The analysis set naming convention is defined by UCSC genome browser. They are:
    - Chromosome are numerically sorted with `chr` prefix. Specific chromosomes are named `chrX` `chrY` and `chrM`.
    - Unlocalized scaffolds are named using `chr[CHROMOSOME]_[ACCESSION]v[VERSION]_random`, where `[ACCESSION]` satisfies GenBank naming conventions.
    - Unplaced scaffolds are named using `chrUn_[ACCESSION]v[VERSION]`.
    - Alternative loci are named using `chr[CHROMOSOME]_[ACCESSION]v[VERSION]_alt`.
    - Fix patches are named using `chr[CHROMOSOME]_[ACCESSION]v[VERSION]_fix`.
    - EBV sequences are placed under `chrEBV`
    - Other decoys named `chrUn_[ACCESSION]v[VERSION]_decoy`.
    - Human Leukocyte Antigen (HLA) sequences are places under `HLA-[HLA_NAME]`.

    See also: {cite}`NCBI_analysis_set` and {cite}`1kg_2015`.

Ensebml Naming Convention
: The naming convention used by Ensembl Genome Browser. They are:
    - Chromosome names without `chr` prefix: `1`, `2` ... `22`, `X`, `Y`, `MT`.
    - Other regions shown by NCBI RefSeq Accession.

GenCode Accesion
: The naming convention used by GenCode project. Is Ensembl convention with `chr` prefix in chromosome names.

## FASTA Samples

Summary:

- AS: Analysis Set
- PA: Primiary Assembly
- TL: Toplevel Assembly
- SM: Soft Mask
- HM/RM: Hard Mask

```{csv-table}
:file: feature.tsv
:delim: '	'
```

```{toctree}
:caption: 'Detailed Description:'
:glob:
:maxdepth: 2

*
```
