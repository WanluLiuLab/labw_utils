# GENCODE Flavor

Description from <https://www.gencodegenes.org/pages/gencode.html>:

> The National Human Genome Research Institute (NHGRI) launched a public research consortium named [ENCODE](http://www.genome.gov/10005107), the Encyclopedia Of DNA Elements, in September 2003, to carry out a project to identify all functional elements in the human genome sequence. After a successful pilot phase on 1% of the genome, the scale-up to the entire genome is now underway. The Wellcome Sanger Institute was [awarded a grant](http://www.sanger.ac.uk/news/view/2007-10-09-researchers-expand-efforts-to-explore-functional-landscape-of-the-human-genome) to carry out a scale-up of the GENCODE project for integrated annotation of gene features.
>
> Having been involved in successfully delivering the definitive annotation of functional elements in the human genome, the GENCODE group were [awarded a second grant](https://www.gencodegenes.org/pages/phase2_goals.html) in 2013 in order to continue their human genome annotation work and expand GENCODE to include annotation of the mouse genome. A [third grant](http://grantome.com/grant/NIH/U41-HG007234-05) was awarded in 2017 for the continued improvement of the annotation of the human and mouse genomes.
>
> The GENCODE gene sets are used by the entire ENCODE consortium and by many other projects (eg. Genotype-Tissue Expression (GTEx), The Cancer Genome Atlas (TCGA), International Cancer Genome Consortium (ICGC), NIH Roadmap Epigenomics Mapping Consortium, Blueprint Epigenome Project, Exome Aggregation Consortium (EXAC), Genome Aggregation Database (gnomAD), 1000 Genomes Project and the Human Cell Atlas (HCA)) as reference gene sets.
>
> The aims of the current GENCODE phase running from 2017 to 2021 are:
>
> - To continue to improve the coverage and accuracy of the GENCODE human and mouse gene sets by enhancing and extending the annotation of all evidence-based gene features in the human genome at a high accuracy, including protein-coding loci with alternatively splices variants, non-coding loci and pseudogenes.
>
> The process to create this annotation involves manual curation, computational analysis and targeted experimental approaches.
>
> The human and mouse GENCODE resources will continue to be available to the research community with regular releases of Ensembl genome browser and the UCSC genome browser will continue to present the current release of the GENCODE gene set.

- Chromosomes: Sorted on numerical flavor, `X`, `Y`, `M`, with `chr` prefix. Other contigs sorted on unknown flavor.
- Coordinate: Related to text order: Incremental on positive strand, decremental on negative strand.
- Quoted Attribute: String quoted.
- `exon_number`: Related to text order: Incremental on positive strand and negative strand.
- `gene_id`: Gene in Ensembl Format
- `gene_name`: Gene in HGNC Format
- `transcript_id`: Transcript in Ensembl Format
- `transcript_name`: Transcript in HGNC Format
- `exon_id`: Exon in Ensembl Format

- Unknown whether there are duplicated keys in attributes.

## Example

Downloaded from <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz>

### List of all chromosomes

`chr1` `chr2` `chr3` `chr4` `chr5` `chr6` `chr7` `chr8` `chr9` `chr10` `chr11` `chr12` `chr13` `chr14` `chr15` `chr16` `chr17` `chr18` `chr19` `chr20` `chr21` `chr22` `chrX` `chrY` `chrM` `GL000009.2` `GL000194.1` `GL000195.1` `GL000205.2` `GL000213.1` `GL000216.2` `GL000218.1` `GL000219.1` `GL000220.1` `GL000225.1` `KQ759759.1` `ML143376.1` `KN538364.1` `ML143355.1` `ML143347.1` `ML143352.1` `KQ759762.1` `ML143375.1` `KQ031383.1` `KN538369.1` `ML143342.1` `ML143350.1` `ML143362.1` `JH159136.1` `ML143357.1` `ML143378.1` `JH159137.1` `ML143365.1` `KZ208924.1` `KQ031387.1` `KV766195.1` `KZ208916.1` `KN538360.1` `KZ208920.1` `KZ208906.1` `KN196484.1` `KN196476.1` `KQ983257.1` `KN196479.1` `KV575245.1` `KZ208917.1` `KZ208911.1` `KN196473.1` `KZ559108.1` `KN196487.1` `KQ759760.1` `KN196475.1` `KV880766.1` `KQ090016.1` `ML143374.1` `KV880764.1` `KN538361.1` `KN196474.1` `ML143360.1` `KZ559109.1` `ML143359.1` `KQ090022.1` `KV766194.1` `KN196478.1` `KZ559104.1` `KN196480.1` `ML143370.1` `KQ090028.1` `KN196481.1` `KN538363.1` `KN538362.1` `KQ031385.1` `KV766192.1` `KQ031388.1` `ML143361.1` `KN538370.1` `KN538373.1` `KZ559113.1` `KV880765.1` `KV766196.1` `KN538371.1` `KQ031384.1` `KN538372.1` `KQ090021.1` `ML143371.1` `KZ559115.1` `KZ208914.1` `KZ208922.1` `ML143373.1` `ML143369.1` `ML143366.1` `ML143367.1` `ML143372.1` `ML143380.1` `ML143377.1` `ML143345.1` `KQ458386.1` `ML143358.1` `KV575244.1` `ML143381.1` `ML143354.1` `ML143344.1` `ML143349.1` `KZ208912.1` `ML143341.1` `KZ208915.1` `KV880768.1` `KN196472.1` `GL383545.1` `GL383546.1` `KI270824.1` `KI270825.1` `KN538368.1` `KZ559111.1` `KI270829.1` `KI270830.1` `KI270831.1` `KI270832.1` `KI270902.1` `KI270903.1` `KZ559110.1` `KI270927.1` `GL877875.1` `GL383550.2` `GL877876.1` `GL383552.1` `KI270904.1` `GL383553.2` `KI270835.1` `GL383551.1` `KI270837.1` `KI270833.1` `KI270834.1` `KZ559112.1` `KI270838.1` `KI270840.1` `KI270842.1` `KI270844.1` `KI270845.1` `KI270846.1` `KI270847.1` `KZ208919.1` `KI270852.1` `KI270848.1` `GL383554.1` `KI270906.1` `GL383555.2` `KI270851.1` `KI270849.1` `KI270905.1` `KI270850.1` `KQ031389.1` `KI270853.1` `GL383556.1` `GL383557.1` `KI270855.1` `KI270856.1` `KQ090027.1` `KQ090026.1` `KZ208921.1` `KI270854.1` `KI270909.1` `KZ559114.1` `GL383563.3` `KI270861.1` `GL383564.2` `GL000258.2` `KI270860.1` `KI270907.1` `KI270862.1` `GL383565.1` `KI270908.1` `KV766198.1` `KI270910.1` `GL383566.1` `JH159146.1` `JH159147.1` `JH159148.1` `KI270857.1` `KI270858.1` `KI270859.1` `KZ559116.1` `GL383567.1` `GL383569.1` `GL383570.1` `GL383571.1` `GL383572.1` `KI270863.1` `KQ458385.1` `KI270912.1` `KI270911.1` `KV575254.1` `KV575246.1` `KV575256.1` `KV575253.1` `KV575252.1` `KV575255.1` `KV575259.1` `KI270917.1` `KI270918.1` `KI270919.1` `KV575247.1` `KV575248.1` `KV575250.1` `KV575249.1` `KV575257.1` `KI270920.1` `KI270921.1` `KI270922.1` `KI270923.1` `KI270929.1` `KI270930.1` `KI270931.1` `KI270932.1` `KI270933.1` `KI270882.1` `KI270883.1` `KI270884.1` `KI270885.1` `KI270886.1` `KI270887.1` `KI270888.1` `KV575258.1` `KV575251.1` `KV575260.1` `KI270889.1` `KI270890.1` `GL000209.2` `KI270891.1` `KI270914.1` `KI270915.1` `KI270916.1` `GL949746.1` `GL949747.2` `GL949748.2` `GL949749.2` `GL949750.2` `GL949751.2` `GL949752.1` `GL949753.2` `GL383573.1` `GL383574.1` `GL383575.2` `KI270866.1` `GL383576.1` `KI270867.1` `KI270865.1` `KI270938.1` `KI270868.1` `KI270760.1` `KI270762.1` `GL383518.1` `KI270759.1` `KI270766.1` `GL383519.1` `KI270761.1` `KQ458382.1` `GL383520.2` `KI270763.1` `KQ458383.1` `KI270765.1` `KQ983255.1` `KQ458384.1` `KV880763.1` `KZ208904.1` `KZ208905.1` `KI270892.1` `GL383577.2` `KI270869.1` `KI270870.1` `KI270871.1` `GL383579.2` `GL383580.2` `GL383581.2` `KI270872.1` `KI270873.1` `KI270874.1` `GL383582.2` `GL383583.2` `KI270875.1` `KI270876.1` `KI270877.1` `KI270878.1` `KI270879.1` `KB663609.1` `KI270928.1` `KN196485.1` `KN196486.1` `KQ458387.1` `KQ458388.1` `KQ759761.1` `KI270769.1` `KI270767.1` `GL383521.1` `KI270772.1` `GL383522.1` `KI270770.1` `KI270893.1` `KI270894.1` `GL582966.2` `KI270773.1` `KI270776.1` `KI270768.1` `KI270774.1` `KQ983256.1` `KZ208907.1` `KZ208908.1` `JH636055.2` `GL383526.1` `KI270779.1` `KI270777.1` `KI270782.1` `KI270783.1` `KI270895.1` `KZ208909.1` `KI270780.1` `KI270924.1` `ML143343.1` `KI270781.1` `KI270934.1` `KZ559105.1` `KI270935.1` `KI270936.1` `KZ559102.1` `KI270937.1` `KZ559103.1` `KI270784.1` `KQ983258.1` `KV766193.1` `GL383527.1` `KI270790.1` `GL383528.1` `GL000257.2` `KI270785.1` `KI270786.1` `KI270788.1` `KI270789.1` `KI270896.1` `KI270925.1` `KQ090014.1` `KQ090015.1` `GL383532.1` `KI270897.1` `GL383531.1` `GL949742.1` `GL339449.2` `KI270795.1` `KI270791.1` `KI270898.1` `KI270792.1` `KI270796.1` `KI270793.1` `KI270794.1` `KN196477.1` `KV575243.1` `GL383533.1` `KB021644.2` `KI270797.1` `KI270798.1` `KI270799.1` `KI270800.1` `KI270801.1` `KI270802.1` `KI270758.1` `GL000250.2` `GL000251.2` `GL000252.2` `GL000253.2` `GL000254.2` `GL000255.2` `GL000256.2` `KI270804.1` `KI270806.1` `GL383534.2` `KI270805.1` `KI270899.1` `KI270809.1` `KI270803.1` `KI270807.1` `KZ559106.1` `KZ208913.1` `KI270808.1` `KI270811.1` `KI270814.1` `KI270810.1` `KI270812.1` `KI270815.1` `KI270813.1` `KI270816.1` `KI270818.1` `KI270817.1` `KI270900.1` `KI270819.1` `KI270901.1` `KI270926.1` `KZ559107.1` `KI270821.1` `KI270822.1` `GL383539.1` `GL383540.1` `GL383541.1` `GL383542.1` `KI270823.1` `KQ090018.1` `KI270880.1` `KI270881.1` `KI270913.1` `KI270442.1` `KI270711.1` `KI270713.1` `KI270721.1` `KI270726.1` `KI270727.1` `KI270728.1` `KI270731.1` `KI270733.1` `KI270734.1` `KI270744.1` `KI270750.1`

### Prefix annotations

```text
##description: evidence-based annotation of the human genome (GRCh38), version 42 (Ensembl 108)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2022-07-20
```

### A transcript on positive strand

```text
chr1	HAVANA	gene	173477330	173488815	.	+	.	gene_id "ENSG00000117592.9"; gene_type "protein_coding"; gene_name "PRDX6"; level 2; hgnc_id "HGNC:16753"; havana_gene "OTTHUMG00000034804.2";
chr1	HAVANA	transcript	173477330	173481971	.	+	.	gene_id "ENSG00000117592.9"; transcript_id "ENST00000460950.1"; gene_type "protein_coding"; gene_name "PRDX6"; transcript_type "processed_transcript"; transcript_name "PRDX6-202"; level 2; transcript_support_level "2"; hgnc_id "HGNC:16753"; havana_gene "OTTHUMG00000034804.2"; havana_transcript "OTTHUMT00000084223.1";
chr1	HAVANA	exon	173477330	173477492	.	+	.	gene_id "ENSG00000117592.9"; transcript_id "ENST00000460950.1"; gene_type "protein_coding"; gene_name "PRDX6"; transcript_type "processed_transcript"; transcript_name "PRDX6-202"; exon_number 1; exon_id "ENSE00001269492.7"; level 2; transcript_support_level "2"; hgnc_id "HGNC:16753"; havana_gene "OTTHUMG00000034804.2"; havana_transcript "OTTHUMT00000084223.1";
chr1	HAVANA	exon	173481326	173481971	.	+	.	gene_id "ENSG00000117592.9"; transcript_id "ENST00000460950.1"; gene_type "protein_coding"; gene_name "PRDX6"; transcript_type "processed_transcript"; transcript_name "PRDX6-202"; exon_number 2; exon_id "ENSE00001880998.1"; level 2; transcript_support_level "2"; hgnc_id "HGNC:16753"; havana_gene "OTTHUMG00000034804.2"; havana_transcript "OTTHUMT00000084223.1";
```

### A transcript on negative strand

```text
chr1	HAVANA	exon	29534	29570	.	-	.	gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_name "WASH7P-201"; exon_number 1; exon_id "ENSE00001890219.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:38034"; ont "PGO:0000005"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1	HAVANA	exon	24738	24891	.	-	.	gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_name "WASH7P-201"; exon_number 2; exon_id "ENSE00003507205.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:38034"; ont "PGO:0000005"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1	HAVANA	exon	18268	18366	.	-	.	gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_name "WASH7P-201"; exon_number 3; exon_id "ENSE00003477500.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:38034"; ont "PGO:0000005"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1	HAVANA	exon	17915	18061	.	-	.	gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_name "WASH7P-201"; exon_number 4; exon_id "ENSE00003565697.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:38034"; ont "PGO:0000005"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1	HAVANA	exon	17606	17742	.	-	.	gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_name "WASH7P-201"; exon_number 5; exon_id "ENSE00003475637.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:38034"; ont "PGO:0000005"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1	HAVANA	exon	17233	17368	.	-	.	gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_name "WASH7P-201"; exon_number 6; exon_id "ENSE00003502542.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:38034"; ont "PGO:0000005"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1	HAVANA	exon	16858	17055	.	-	.	gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_name "WASH7P-201"; exon_number 7; exon_id "ENSE00003553898.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:38034"; ont "PGO:0000005"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1	HAVANA	exon	16607	16765	.	-	.	gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_name "WASH7P-201"; exon_number 8; exon_id "ENSE00003621279.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:38034"; ont "PGO:0000005"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1	HAVANA	exon	15796	15947	.	-	.	gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_name "WASH7P-201"; exon_number 9; exon_id "ENSE00002030414.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:38034"; ont "PGO:0000005"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1	HAVANA	exon	15005	15038	.	-	.	gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_name "WASH7P-201"; exon_number 10; exon_id "ENSE00001935574.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:38034"; ont "PGO:0000005"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
chr1	HAVANA	exon	14404	14501	.	-	.	gene_id "ENSG00000227232.5"; transcript_id "ENST00000488147.1"; gene_type "unprocessed_pseudogene"; gene_name "WASH7P"; transcript_type "unprocessed_pseudogene"; transcript_name "WASH7P-201"; exon_number 11; exon_id "ENSE00001843071.1"; level 2; transcript_support_level "NA"; hgnc_id "HGNC:38034"; ont "PGO:0000005"; tag "basic"; tag "Ensembl_canonical"; havana_gene "OTTHUMG00000000958.1"; havana_transcript "OTTHUMT00000002839.1";
```

## Official Readme: Format description of GENCODE GTF

Downloaded from <https://www.gencodegenes.org/pages/data_format.html>, <https://www.gencodegenes.org/pages/tags.html>

### A. TAB-separated standard GTF columns

chromosome name
:`chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M}` or GRC accession. Scaffolds, patches and haplotypes names correspond to their GRC accessions. Please note that these are different from the Ensembl names.

annotation source
: `ENSEMBL` or `HAVANA`

feature type
: `gene` `transcript` `exon` `CDS` `UTR` `start_codon` `stop_codon` `Selenocysteine`

genomic start location
: integer-value (1-based)

genomic end location
: integer-value

score
: not used

genomic strand
: `+` `-`

genomic phase (for CDS features)
: `0` `1` `2` `.`

additional information as key-value pairs (see below)

### B. Key-value pairs in 9th column (format: key "value"; )

% To be refactored

#### B.1. Mandatory fields

<table class="striped"
data-summary="Mandatory fields and their format in the 9th column of GTF">
<colgroup>
<col style="width: 25%" />
<col style="width: 25%" />
<col style="width: 25%" />
<col style="width: 25%" />
</colgroup>
<thead>
<tr class="header">
<th>key name</th>
<th>feature type(s)</th>
<th>value format</th>
<th>release</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>gene_id</td>
<td>all</td>
<td>ENSGXXXXXXXXXXX.X <sup>b,c</sup> _X<sup>g</sup></td>
<td>all</td>
</tr>
<tr class="even tint">
<td>transcript_id <sup>d</sup></td>
<td>all except gene</td>
<td>ENSTXXXXXXXXXXX.X <sup>b,c</sup> _X<sup>g</sup></td>
<td>all</td>
</tr>
<tr class="odd">
<td>gene_type</td>
<td>all</td>
<td><a href="/pages/biotypes.html">list of biotypes</a></td>
<td>all</td>
</tr>
<tr class="even tint">
<td>gene_status <sup>e</sup></td>
<td>all</td>
<td>{KNOWN, NOVEL, PUTATIVE}</td>
<td>until 25 and M11</td>
</tr>
<tr class="odd">
<td>gene_name</td>
<td>all</td>
<td>string</td>
<td>all</td>
</tr>
<tr class="even tint">
<td>transcript_type <sup>d</sup></td>
<td>all except gene</td>
<td><a href="/pages/biotypes.html">list of biotypes</a></td>
<td>all</td>
</tr>
<tr class="odd">
<td>transcript_status<sup>d,e</sup></td>
<td>all except gene</td>
<td>{KNOWN, NOVEL, PUTATIVE}</td>
<td>until 25 and M11</td>
</tr>
<tr class="even tint">
<td>transcript_name <sup>d</sup></td>
<td>all except gene</td>
<td>string</td>
<td>all</td>
</tr>
<tr class="odd">
<td>exon_number <sup>f</sup></td>
<td>all except gene/transcript/Selenocysteine</td>
<td>integer (exon position in the transcript from its 5' end)</td>
<td>all</td>
</tr>
<tr class="even tint">
<td>exon_id <sup>f</sup></td>
<td>all except gene/transcript/Selenocysteine</td>
<td>ENSEXXXXXXXXXXX.X <sup>b</sup> _X<sup>g</sup></td>
<td>all</td>
</tr>
<tr class="odd">
<td>level</td>
<td>all</td>
<td>1 (verified loci),<br />
2 (manually annotated loci),<br />
3 (automatically annotated loci)</td>
<td>all</td>
</tr>
</tbody>
</table>

<sup>b</sup> From version 7 the gene/transcript version number was
appended to gene and transcript ids (eg. ENSG00000160087.16).

<sup>c</sup> Gene and trancript ids on the chrY PAR regions have
"\_PAR_Y" appended (from release 25), or are in the format
ENSGRXXXXXXXXXX and ENSTRXXXXXXXXXX (until release 24) to avoid
redundancy.

<sup>d</sup> Until releases 21 and M4, the gene lines included
transcript attributes.

<sup>e</sup> The 'gene_status' and 'transcript_status' attributes were
removed after releases 25 (human) and M11 (mouse).

<sup>f</sup> Except in gene and transcript lines.

<sup>g</sup> In the annotation mapped back to GRCh37, mapping versions
are appended to the identifiers (eg. ENSG00000228327.3_2).

#### B.2. Optional fields

<table class="striped"
data-summary="Optional fields and their format in the 9th column of GTF">
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<thead>
<tr class="header">
<th>key name</th>
<th>value format</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>tag</td>
<td>part of a special set [*]: <a href="/pages/tags.html">list of
tags</a></td>
</tr>
<tr class="even tint">
<td>ccdsid</td>
<td>official CCDS id [*]; CCDS*</td>
</tr>
<tr class="odd">
<td>havana_gene</td>
<td>gene id in the havana db [0,1]; OTTHUMGXXXXXXXXXXX.X</td>
</tr>
<tr class="even tint">
<td>havana_transcript</td>
<td>transcript id in the havana db [0,1] ; OTTHUMTXXXXXXXXXXX.X</td>
</tr>
<tr class="odd">
<td>protein_id</td>
<td>ENSPXXXXXXXXXXX.X [0,1]</td>
</tr>
<tr class="even tint">
<td>ont</td>
<td>pseudogene (or other) ontology ids [*]; {PGO:0000004 and
others}</td>
</tr>
<tr class="odd">
<td>transcript_support_level</td>
<td>{1,2,3,4,5,NA} [0,1]<br />
transcripts are scored according to how well mRNA and EST alignments
match over its full length:<br />
1 (all splice junctions of the transcript are supported by at least one
non-suspect mRNA),<br />
2 (the best supporting mRNA is flagged as suspect or the support is from
multiple ESTs),<br />
3 (the only support is from a single EST),<br />
4 (the best supporting EST is flagged as suspect),<br />
5 (no single transcript supports the model structure),<br />
NA (the transcript was not analyzed)</td>
</tr>
<tr class="even tint">
<td>remap_status<br />
remap_original_id<br />
remap_original_location<br />
remap_num_mappings<br />
remap_target_status<br />
remap_substituted_missing_target</td>
<td><a href="/human/grch37_mapped_releases.html#attrib">Mapping
attributes</a> [0,1] - only for GRCh38 annotation lifted back to
GRCh37.</td>
</tr>
<tr class="odd">
<td>hgnc_id</td>
<td>HGNC id in human [0,1]; HGNC:*</td>
</tr>
<tr class="even tint">
<td>mgi_id</td>
<td>MGI id in mouse [0,1]; MGI:*</td>
</tr>
</tbody>
</table>

Number of occurrences: \[\*\] - zero or multiple, \[0,1\] - zero or one

### Example GTF lines

(Omitted)

### Examples for fetching specific parts from the file \[Unix command line\]

1. Get all "gene" lines:

    ```shell
    awk '{if($3=="gene"){print $0}}' gencode.gtf
    ```

2. Get all "protein-coding transcript" lines:

    ```shell
    awk '{if($3=="transcript" && $20=="\"protein_coding\";"){print $0}}' gencode.gtf
    ```

3. Get level 1 & 2 annotation (manually annotated) only:

    ```shell
    awk '{if($0~"level (1|2);"){print $0}}' gencode.gtf
    ```

### Example for parsing the file \[Perl\]

```perl
#!/usr/bin/perl

use strict;

my $gencode_file = "gencode.v23.annotation.gtf";
open(IN, "<$gencode_file") or die "Can't open $gencode_file.\n";
my %all_genes;
while(<IN>){
    next if(/^##/); #ignore header
    chomp;
    my %attribs = ();
    my ($chr, $source, $type, $start, $end, $score,
        $strand, $phase, $attributes) = split("\t");
    #store nine columns in hash
    my %fields = (
        chr        => $chr,
        source     => $source,
        type       => $type,
        start      => $start,
        end        => $end,
        score      => $score,
        strand     => $strand,
        phase      => $phase,
        attributes => $attributes,
    );
    my @add_attributes = split(";", $attributes);
    # store ids and additional information in second hash
    foreach my $attr ( @add_attributes ) {
        next unless $attr =~ /^\s*(.+)\s(.+)$/;
        my $c_type  = $1;
        my $c_value = $2;
        $c_value =~ s/\"//g;
        if($c_type  && $c_value){
            if(!exists($attribs{$c_type})){
                $attribs{$c_type} = [];
            }
            push(@{ $attribs{$c_type} }, $c_value);
        }
    }
    #work with the information from the two hashes...
    #eg. store them in a hash of arrays by gene_id:
    if(!exists($all_genes{$attribs{'gene_id'}->[0]})){
        $all_genes{$attribs{'gene_id'}->[0]} = [];
    }
    push(@{ $all_genes{$attribs{'gene_id'}->[0]} }, \%fields);
}
print "Example entry ENSG00000183186.7: ".
    $all_genes{"ENSG00000183186.7"}->[0]->{"type"}.", ".
    $all_genes{"ENSG00000183186.7"}->[0]->{"chr"}." ".
    $all_genes{"ENSG00000183186.7"}->[0]->{"start"}."-".
    $all_genes{"ENSG00000183186.7"}->[0]->{"end"}."\n";
```

### GTF Parsers in Other Programming Language

A number of programming languages already have GTF parsers developed by
third party libraries. We have listed a number of these below and should
be used in preference to writing your own parser.

Java: [BioJava](https://biojava.org/)
([JavaDoc](https://biojava.org/docs/api5.0.0/org/biojava/nbio/genome/parsers/gff/GFF3Reader.html))

Perl: [BioPerl](https://bioperl.org/)
([metacpan](https://metacpan.org/pod/Bio::FeatureIO::gtf))

Python: [gtfparse](https://github.com/openvax/gtfparse),
[gffutils](https://github.com/daler/gffutils)

Ruby: [germ](https://github.com/mountetna/germ)

### Tag description in GENCODE

The following tags can be found in the GENCODE GTF/GFF3 files. Read more about the [GTF file format](https://www.gencodegenes.org/pages/data_format.html)

3_nested_supported_extension
: 3' end extended based on RNA-seq data.

3_standard_supported_extension
: 3' end extended based on RNA-seq data.

454_RNA_Seq_supported
: annotated based on RNA-seq data.

5_nested_supported_extension
: 5' end extended based on RNA-seq data.

5_standard_supported_extension
: 5' end extended based on RNA-seq data.

alternative_3\_UTR
: shares an identical CDS but has alternative 5' UTR with respect to a reference variant.

alternative_5\_UTR
: shares an identical CDS but has alternative 3' UTR with respect to a reference variant.

appris_principal_1
: (This flag corresponds to the older flag "appris_principal") Where the transcript expected to code for the main functional isoform based solely on the core modules in the APPRIS database. The APPRIS core modules map protein structural and functional information and cross-species conservation to the annotated variants.

appris_principal_2
: (This flag corresponds to the older flag "appris_candidate_ccds") Where the APPRIS core modules are unable to choose a clear principal variant (approximately 25% of human protein coding genes), the database chooses two or more of the CDS variants as "candidates" to be the principal variant. If one (but no more than one) of these candidates has a distinct CCDS identifier it is selected as the principal variant for that gene. A CCDS identifier shows that there is consensus between RefSeq and GENCODE/Ensembl for that variant, guaranteeing that the variant has cDNA support.

appris_principal_3
: Where the APPRIS core modules are unable to choose a clear principal variant and there more than one of the variants have distinct CCDS identifiers, APPRIS selects the variant with lowest CCDS identifier as the principal variant. The lower the CCDS identifier, the earlier it was annotated. Consensus CDS annotated earlier are likely to have more cDNA evidence. Consecutive CCDS identifiers are not included in this flag, since they will have been annotated in the same release of CCDS. These are distinguished with the next flag.

appris_principal_4
: (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_ccds") Where the APPRIS core modules are unable to choose a clear principal CDS and there is more than one variant with a distinct (but consecutive) CCDS identifiers, APPRIS selects the longest CCDS isoform as the principal variant.

appris_principal_5
: (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_seq") Where the APPRIS core modules are unable to choose a clear principal variant and none of the candidate variants are annotated by CCDS, APPRIS selects the longest of the candidate isoforms as the principal variant.

appris_alternative_1
: Candidate transcript(s) models that are conserved in at least three tested non-primate species.

appris_alternative_2
: Candidate transcript(s) models that appear to be conserved in fewer than three tested non-primate species.

appris_principal
: transcript expected to code for the main functional isoform based on a range of protein features (APPRIS pipeline).

appris_candidate
: where there is no single 'appris_principal' variant the main functional isoform will be translated from one of the 'appris_candidate' genes.

appris_candidate_ccds
: the "appris_candidate" transcript that has an unique CCDS.

appris_candidate_highest_score
: where there is no 'appris_principal' variant, the candidate with highest APPRIS score is selected as the primary variant.

appris_candidate_longest
: where there is no 'appris_principal' variant, the longest of the 'appris_candidate' variants is selected as the primary variant.

appris_candidate_longest_ccds
: the "appris_candidate" transcripts where there are several CCDS, in this case APPRIS labels the longest CCDS.

appris_candidate_longest_seq
: where there is no "appris_candidate_ccds" or "appris_candidate_longest_ccds" variant, the longest protein of the "appris_candidate" variants is selected as the primary variant.

artifactual_duplication
: annotated on an artifactual duplicate region of the genome assembly.

basic
: identifies a subset of representative transcripts for each gene; prioritises full-length protein coding transcripts over partial or non-protein coding transcripts within the same gene, and intends to highlight those transcripts that will be useful to the majority of users.

bicistronic
: transcript contains two confidently annotated CDSs. Support may come from eg proteomic data, cross-species conservation or published experimental work.

CAGE_supported_TSS
: transcript 5' end overlaps ENCODE or Fantom CAGE cluster.

CCDS
: member of the consensus CDS gene set, confirming coding regions between ENSEMBL, UCSC, NCBI and HAVANA.

cds_end_NF
: the coding region end could not be confirmed.

cds_start_NF
: the coding region start could not be confirmed.

dotter_confirmed
: transcript QC checked using dotplot to identify features eg splice junctions, end of homology.

downstream_ATG
: an upstream ATG is used where a downstream ATG seems more evolutionary conserved.

Ensembl_canonical
: most representative transcript of the gene. This will be the MANE_Select transcript if there is one, or a transcript chosen by an Ensembl algorithm otherwise.

exp_conf
: transcript was tested and confirmed experimentally.

fragmented_locus
: locus consists of non-overlapping transcript fragments either because of genome assembly issues (i.e., gaps or mis-assemblies), or because supporting transcripts (e.g., from another species) cannot be completely mapped, or because the supporting transcripts are non-overlapping end pairs (i.e., 5' and 3' ESTs from a single cDNA).

inferred_exon_combination
: transcript model contains all possible in-frame exons supported by homology, experimental evidence or conservation, but the exon combination is not directly supported by a single piece of evidence and may not be biological. Used for large genes with repetitive exons (e.g. titin (TTN)) to represent all the exons individual transcript variants can pool from.

inferred_transcript_model
: transcript model is not supported by a single piece of transcript evidence. May be supported by multiple fragments of transcript evidence or by combining different evidence sources e.g. protein homology, RNA-seq data, published experimental data.

low_sequence_quality
: transcript supported by transcript evidence that, while ampping best-in-genome, shows regions of poor sequence quality.

mRNA_end_NF
: the mRNA end could not be confirmed.

mRNA_start_NF
: the mRNA start could not be confirmed.

MANE_Select
: the transcript belongs to the MANE Select data set. The Matched Annotation from NCBI and EMBL-EBI project (MANE) is a collaboration between Ensembl-GENCODE and RefSeq to select a default transcript per human protein coding locus that is representative of biology, well-supported, expressed and conserved. This transcript set matches GRCh38 and is 100% identical between RefSeq and Ensembl-GENCODE for 5' UTR, CDS, splicing and 3' UTR.

MANE_Plus_Clinical
: the transcript belongs to the MANE Plus Clinical data set. Within the MANE project, these are additional transcripts per locus necessary to support clinical variant reporting, for example transcripts containing known pathogenic or likely pathogenic clinical variants not reportable using the MANE Select data set. This transcript set matches GRCh38 and is 100% identical between RefSeq and Ensembl-GENCODE for 5' UTR, CDS, splicing and 3' UTR.

NAGNAG_splice_site
: in-frame type of variation where, at the acceptor site, some variants splice after the first AG and others after the second AG.

ncRNA_host
: the locus is a host for small non-coding RNAs.

nested_454_RNA_Seq_supported
: annotated based on RNA-seq data.

NMD_exception
: the transcript looks like it is subject to NMD but publications, experiments or conservation support the translation of the CDS.

NMD_likely_if_extended
: codon if the transcript were longer but cannot currently be annotated as NMD as does not fulfil all criteria - most commonly lack of an intron downstream of the stop codon.

non_ATG_start
: the CDS has a non-ATG start and its validity is supported by publication or conservation.

non_canonical_conserved
: the transcript has a non-canonical splice site conserved in other species.

non_canonical_genome_sequence_error
: the transcript has a non-canonical splice site explained by a genomic sequencing error.

non_canonical_other
: the transcript has a non-canonical splice site explained by other reasons.

non_canonical_polymorphism
: the transcript has a non-canonical splice site explained by a SNP.

non_canonical_TEC
: the transcript has a non-canonical splice site that needs experimental confirmation.

non_canonical_U12
: the transcript has a non-canonical splice site explained by a U12 intron (i.e. AT-AC splice site).

non_submitted_evidence
: a splice variant for which supporting evidence has not been submitted to databases, i.e. the model is based on literature or collaborator evidence.

not_best_in_genome_evidence
: a transcript is supported by evidence from same species paralogous loci.

not_organism_supported
: evidence from other species was used to build model.

orphan
: protein-coding locus with no paralogues or orthologs.

overlapping_locus
: exon(s) of the locus overlap exon(s) of a readthrough transcript or a transcript belonging to another locus.

overlapping_uORF
: a low confidence upstream ATG existing in other coding variant would lead to NMD in this trancript, that uses the high confidence canonical downstream ATG.

PAR
: annotation in the pseudo-autosomal region, which is duplicated between chromosomes X and Y.

pseudo_consens
: member of the pseudogene set predicted by YALE, UCSC and HAVANA.

readthrough_gene
: protein-coding gene that has a readthrough transcript.

readthrough_transcript
: a transcript that overlaps two or more independent loci but is considered to belong to a third, separate locus.

reference_genome_error
: locus overlaps a sequence error or an assembly error in the reference genome that affects its annotation (e.g., 1 or 2bp insertion/deletion, substitution causing premature stop codon). The main effect is that affected transcripts that would have had a CDS are currently annotated without one.

retained_intron_CDS
: internal intron of CDS portion of transcript is retained.

retained_intron_final
: final intron of CDS portion of transcript is retained.

retained_intron_first
: first intron of CDS portion of transcript is retained.

retrogene
: protein-coding locus created via retrotransposition.

RNA_Seq_supported_only
: transcript supported by RNAseq data and not supported by mRNA or EST evidence.

RNA_Seq_supported_partial
: transcript annotated based on mixture of RNA-seq data and EST/mRNA/protein evidence.

RP_supported_TIS
: transcript that contains a CDS that has a translation initiation site supported by Ribosomal Profiling data.

seleno
: contains a selenocysteine.

semi_processed
: a processed pseudogene with one or more introns still present. These are likely formed through the retrotransposition of a retained intron transcript.

sequence_error
: transcript contains at least 1 non-canonical splice junction that is associated with a known or novel genome sequence error.

stop_codon_readthrough
: Transcript whose coding sequence contains an internal stop codon that does not cause the translation termination

TAGENE
: Transcript created or extended using assembled RNA-seq long reads.

upstream_ATG
: an upstream ATG exists when a downstream ATG is better supported.

upstream_uORF
: a low confidence upstream ATG existing in other coding variant would lead to NMD in this trancript, that uses the high confidence canonical downstream ATG.

(C) 2022 [EMBL-EBI](http:/www.ebi.ac.uk).
