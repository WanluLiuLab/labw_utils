# On Standardization \& Sorting of GTF

## Introduction

Gene Transfer Format (GFF) is a commonly-used format for genomic features like genes, transcripts and exons. However, there are differences in formats, record orders, naming conventions of GTF produced by different genomic centers. This document clearifies basic concepts, describes the differences with help of examples provided by corresponding genomic centers.

Before discussing the differences, we would like to discuss the definitions of terms used in this documentation. They are:

Genes
: A biological meaningful region which consists of several isoforms (represented as transcripts on GTF).

Transcripts
: A mRNA molecule made from the genomic DNA. It should contain one or more exons.

### Gene Naming Conventions

There are multiple ways in representing a gene, a transcript or a protein. For example, HGNC Gene Symbol [A1BG-AS1](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:37133) and Ensembl Gene [ENSG00000268895](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000268895) is the same gene. Commonly seen naming conventions are as follows:

- [HUGO Gene Nomenclature Committee (HGNC)](https://www.genenames.org/), which is a group under [Human Genome Organization (HUGO)](https://www.hugo-international.org/), is responsible for approving unique symbols and names for human loci, including protein coding genes, ncRNA genes and pseudogenes, to allow unambiguous scientific communication. It defines common names (e.g., [BRCA1 DNA repair associated](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:1100)), id (e.g., [HGNC:1100](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:1100)) and symbols (e.g., [BRCA1](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:1100)) for a gene for human and mice.
- [Ensembl](https://www.ensembl.org) is a genome browser for vertebrate genomes that supports research in comparative genomics, evolution, sequence variation and transcriptional regulation. It provides accession for genes (e.g., [ENSG00000012048](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000012048)), transcripts (e.g., [ENSG00000012048](https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000012048)), exons (e.g., ENSE00001814242) and proteins (e.g., [ENSG00000012048](https://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;g=ENSG00000012048)).
- [Universal Protein Resource (UniProt)](https://www.uniprot.org), which is a comprehensive resource for protein sequence and annotation data, provides accession number for proteins (e.g., [P38398](https://www.uniprot.org/uniprot/P38398)).
- The [NCBI RefSeq database](https://www.ncbi.nlm.nih.gov/refseq/) from [National Center for Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/) of United States of America produces ID to Genes (e.g., [672](https://www.ncbi.nlm.nih.gov/gene/672)), Transcripts (e.g., [NM_007294](https://www.ncbi.nlm.nih.gov/nuccore/NM_007294), or [NR_024540](https://www.ncbi.nlm.nih.gov/nuccore/NR_024540)), Coding Sequence (CDS) (e.g., [HSU14680](https://www.ncbi.nlm.nih.gov/nuccore/U14680)), Consensus CDS (CCDS) (e.g., [CCDS11453.1](https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=CCDS11453.1)) and Proteins (e.g., [6552299](https://www.ncbi.nlm.nih.gov/protein/6552299)).
- [UCSC Genome Browser](http://www.genome.ucsc.edu/) by University of California, Santa Cruz is another Ensembl-like genome browser. Alongside providing naming conversions, it also provides accession to genes (e.g., [uc002ict.4](http://genome.cse.ucsc.edu/cgi-bin/hgGene?org=Human&hgg_chrom=none&hgg_type=knownGene&hgg_gene=uc002ict.4))

Other formats may include:

- [Mouse Genome Informatics (MGI)](http://www.informatics.jax.org/), which the international database resource for the laboratory mouse, provides accession for genes (e.g., [MGI:104537](http://www.informatics.jax.org/marker/MGI:104537)).
- [Rat Genome Database (RGD)](https://rgd.mcw.edu/)
- [Online Mendelian Inheritance in Man (OMIM)](https://omim.org/), which is a comprehensive, authoritative compendium of human genes and genetic phenotypes, provides accesion for genes (e.g., [113705](https://www.omim.org/entry/113705)).
- [WormBase](https://wormbase.org)
- [FlyBase](https://flybase.org)

## External sorting algorithms

`bedtools`: Sorted by chromosome (alphabetical or self-provided) and coordinate (regardless of strand). Exon number unchanged.

`gffread`: Sorted by chromosome (alphabetical or self-provided) and coordinate (regardless of strand). Exon number depleted.

## GTF Format Specifications

```{toctree}
:caption: 'Contents:'
:glob:
:maxdepth: 2

*
```
