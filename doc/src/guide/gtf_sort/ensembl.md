# Ensembl Flavor

- Chromosomes: Sorted on unknown flavor, does not have `chr` prefix.
- Coordinate: Related to text order: Incremental on positive strand, decremental on negative strand.
- Quoted Attribute: All quoted.
- `exon_number`: Related to text order: Incremental on positive strand and negative strand.
- `gene_id`: Gene in Ensembl Format
- `gene_name`: Gene in HGNC Format
- `transcript_id`: Transcript in Ensembl Format
- `transcript_name`: Transcript in HGNC Format
- `exon_id`: Exon in Ensembl Format

- Duplicated keys in attributes. See readme attached below.

## Example

Downloaded from: <https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz>

### List of all chromosomes

`1` `2` `3` `4` `5` `6` `7` `X` `8` `9` `11` `10` `12` `13` `14` `15` `16` `17` `18` `20` `19` `Y` `22` `21` `MT` `KI270728.1` `KI270727.1` `KI270442.1` `GL000225.1` `GL000009.2` `GL000194.1` `GL000205.2` `GL000195.1` `KI270733.1` `GL000219.1` `GL000216.2` `KI270744.1` `KI270734.1` `GL000213.1`

### Prefix annotations

```text
#!genome-build GRCh38.p13
#!genome-version GRCh38
#!genome-date 2013-12
#!genome-build-accession GCA_000001405.28
#!genebuild-last-updated 2021-03
```

### A transcript on positive strand

```text
1	ensembl_havana	gene	173477330	173488815	.	+	.	gene_id "ENSG00000117592"; gene_version "9"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
1	ensembl_havana	transcript	173477335	173488815	.	+	.	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	exon	173477335	173477492	.	+	.	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "1"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; exon_id "ENSE00001887950"; exon_version "2"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	CDS	173477398	173477492	.	+	0	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "1"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; protein_id "ENSP00000342026"; protein_version "5"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	start_codon	173477398	173477400	.	+	0	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "1"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	exon	173481326	173481482	.	+	.	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "2"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; exon_id "ENSE00000814648"; exon_version "1"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	CDS	173481326	173481482	.	+	1	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "2"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; protein_id "ENSP00000342026"; protein_version "5"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	exon	173485361	173485507	.	+	.	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "3"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; exon_id "ENSE00003615087"; exon_version "1"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	CDS	173485361	173485507	.	+	0	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "3"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; protein_id "ENSP00000342026"; protein_version "5"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	exon	173486255	173486401	.	+	.	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "4"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; exon_id "ENSE00003682006"; exon_version "1"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	CDS	173486255	173486401	.	+	0	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "4"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; protein_id "ENSP00000342026"; protein_version "5"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	exon	173487735	173488815	.	+	.	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "5"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; exon_id "ENSE00001161142"; exon_version "7"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	CDS	173487735	173487860	.	+	0	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "5"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; protein_id "ENSP00000342026"; protein_version "5"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	stop_codon	173487861	173487863	.	+	0	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; exon_number "5"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	five_prime_utr	173477335	173477397	.	+	.	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
1	ensembl_havana	three_prime_utr	173487864	173488815	.	+	.	gene_id "ENSG00000117592"; gene_version "9"; transcript_id "ENST00000340385"; transcript_version "6"; gene_name "PRDX6"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRDX6-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS1307"; tag "basic"; transcript_support_level "1 (assigned to previous version 5)";
```

### A transcript on negative strand

```text
1	havana	exon	29534	29570	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000488147"; transcript_version "1"; exon_number "1"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00001890219"; exon_version "1"; tag "basic"; transcript_support_level "NA";
1	havana	exon	24738	24891	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000488147"; transcript_version "1"; exon_number "2"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00003507205"; exon_version "1"; tag "basic"; transcript_support_level "NA";
1	havana	exon	18268	18366	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000488147"; transcript_version "1"; exon_number "3"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00003477500"; exon_version "1"; tag "basic"; transcript_support_level "NA";
1	havana	exon	17915	18061	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000488147"; transcript_version "1"; exon_number "4"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00003565697"; exon_version "1"; tag "basic"; transcript_support_level "NA";
1	havana	exon	17606	17742	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000488147"; transcript_version "1"; exon_number "5"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00003475637"; exon_version "1"; tag "basic"; transcript_support_level "NA";
1	havana	exon	17233	17368	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000488147"; transcript_version "1"; exon_number "6"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00003502542"; exon_version "1"; tag "basic"; transcript_support_level "NA";
1	havana	exon	16858	17055	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000488147"; transcript_version "1"; exon_number "7"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00003553898"; exon_version "1"; tag "basic"; transcript_support_level "NA";
1	havana	exon	16607	16765	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000488147"; transcript_version "1"; exon_number "8"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00003621279"; exon_version "1"; tag "basic"; transcript_support_level "NA";
1	havana	exon	15796	15947	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000488147"; transcript_version "1"; exon_number "9"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00002030414"; exon_version "1"; tag "basic"; transcript_support_level "NA";
1	havana	exon	15005	15038	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000488147"; transcript_version "1"; exon_number "10"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00001935574"; exon_version "1"; tag "basic"; transcript_support_level "NA";
1	havana	exon	14404	14501	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000488147"; transcript_version "1"; exon_number "11"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00001843071"; exon_version "1"; tag "basic"; transcript_support_level "NA";
```

## Official Readme

Downloaded from <https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/README>

The GTF (General Transfer Format) is an extension of GFF version 2 and used to represent transcription models. GFF (General Feature Format) consists of one line per feature, each containing 9 columns of data. 

### Fields

Fields are tab-separated. Also, all but the final field in each 
feature line must contain a value; "empty" columns are denoted 
with a '.'

- `seqname`: name of the chromosome or scaffold; chromosome names without a `chr`
- `source`: name of the program that generated this feature, or the data source (database or project name)
- `feature`: feature type name. Current allowed features are {gene, transcript, exon, CDS, Selenocysteine, start_codon, stop_codon and UTR}
- `start`: start position of the feature, with sequence numbering starting at 1.
- `end`: end position of the feature, with sequence numbering starting at 1.
- `score`: a floating point value indiciating the score of a feature
- `strand`: defined as `+` (forward) or `-` (reverse).
- `frame`:  one of `0`, `1` or `2`. Frame indicates the number of base pairs before you encounter a full codon. All values are given with relation to the 5' end.
  - `0` indicates the feature begins with a whole codon.
  - `1` indicates there is an extra base (the 3rd base of the prior codon) at the start of this feature.
  - `2` indicates there are two extra bases (2nd and 3rd base of the prior exon) before the first codon.
- `attribute`: a semicolon-separated list of tag-value pairs (separated by a space), providing additional information about each feature. A key can be repeated multiple times.

### Attributes

The following attributes are available. All attributes are semi-colon
separated pairs of keys and values.

- `gene_id`: The stable identifier for the gene
- `gene_version`: The stable identifier version for the gene
- `gene_name`: The official symbol of this gene
- `gene_source`: The annotation source for this gene
- `gene_biotype`: The biotype of this gene
- `transcript_id`: The stable identifier for this transcript
- `transcript_version`: The stable identifier version for this transcript
- `transcript_name`: The symbold for this transcript derived from the gene name
- `transcript_source`: The annotation source for this transcript
- `transcript_biotype`: The biotype for this transcript
- `exon_id`: The stable identifier for this exon
- `exon_version`: The stable identifier version for this exon
- `exon_number`: Position of this exon in the transcript
- `ccds_id`: CCDS identifier linked to this transcript
- `protein_id`: Stable identifier for this transcript's protein
- `protein_version`: Stable identifier version for this transcript's protein
- `tag`: A collection of additional key value tags
- `transcript_support_level`: Ranking to assess how well a transcript is supported (from 1 to 5)

### Tags

Tags are additional flags used to indicate attibutes of the transcript.

- `CCDS`: Flags this transcript as one linked to a CCDS record
- `seleno`: Flags this transcript has a Selenocysteine edit. Look for the Selenocysteine
feature for the position of this on the genome
- `cds_end_NF`: the coding region end could not be confirmed
- `cds_start_NF`: the coding region start could not be confirmed
- `mRNA_end_NF`: the mRNA end could not be confirmed
- `mRNA_start_NF`: the mRNA start could not be confirmed.
- `basic`: the transcript is part of the gencode basic geneset

### Comments

Lines may be commented out by the addition of a single `#` character at the start. These lines should be ignored by your parser.

### Pragmas/Metadata

GTF files can contain meta-data. In the case of experimental meta-data these are noted by a `#!`. Those which are stable are noted by a `##`. Meta data is a single key, a space and then the value. Current meta data keys are:

- `genome-build`: Build identifier of the assembly e.g. GRCh37.p11
- `genome-version`: Version of this assembly e.g. GRCh37
- `genome-date`: The date of this assembly's release e.g. 2009-02
- `genome-build-accession`: The accession and source of this accession e.g. NCBI:GCA_000001405.14
- `genebuild-last-updated`: The date of the last genebuild update e.g. 2013-09
