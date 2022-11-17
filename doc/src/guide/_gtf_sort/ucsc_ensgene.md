# UCSC Flavor of `ensGene`

- Chromosomes: Sorted on alphabetical order with `chr` prefix, where `chr10` after `chr1`.
- Coordinate: Related to text order: Incremental on positive and negative strand.
- Quoted Attribute: All quoted.
- `exon_number`: Related to text order: Incremental on positive and negative strand.
- `gene_id`: Gene in Ensembl Format
- `gene_name`: Gene in Ensembl Format
- `transcript_id`: Transcript in Ensembl Format
- `transcript_name`: NA
- `exon_id`: Exon in Ensembl Format

## Examples

### List of all chromosomes

`chr1` `chr10` `chr10_GL383545v1_alt` `chr10_GL383546v1_alt` `chr10_KI270825v1_alt` `chr11` `chr11_JH159136v1_alt` `chr11_JH159137v1_alt` `chr11_KI270721v1_random` `chr11_KI270829v1_alt` `chr11_KI270830v1_alt` `chr11_KI270831v1_alt` `chr11_KI270832v1_alt` `chr11_KI270902v1_alt` `chr11_KI270903v1_alt` `chr11_KI270927v1_alt` `chr12` `chr12_GL383550v2_alt` `chr12_GL383551v1_alt` `chr12_GL383552v1_alt` `chr12_GL383553v2_alt` `chr12_GL877875v1_alt` `chr12_GL877876v1_alt` `chr12_KI270833v1_alt` `chr12_KI270834v1_alt` `chr12_KI270835v1_alt` `chr12_KI270837v1_alt` `chr12_KI270904v1_alt` `chr13` `chr13_KI270838v1_alt` `chr13_KI270840v1_alt` `chr14` `chr14_GL000009v2_random` `chr14_GL000194v1_random` `chr14_GL000225v1_random` `chr14_KI270722v1_random` `chr14_KI270723v1_random` `chr14_KI270724v1_random` `chr14_KI270726v1_random` `chr14_KI270844v1_alt` `chr14_KI270845v1_alt` `chr14_KI270846v1_alt` `chr14_KI270847v1_alt` `chr15` `chr15_GL383554v1_alt` `chr15_GL383555v2_alt` `chr15_KI270727v1_random` `chr15_KI270848v1_alt` `chr15_KI270849v1_alt` `chr15_KI270850v1_alt` `chr15_KI270851v1_alt` `chr15_KI270852v1_alt` `chr15_KI270905v1_alt` `chr15_KI270906v1_alt` `chr16` `chr16_GL383556v1_alt` `chr16_GL383557v1_alt` `chr16_KI270728v1_random` `chr16_KI270853v1_alt` `chr16_KI270854v1_alt` `chr16_KI270855v1_alt` `chr16_KI270856v1_alt` `chr17` `chr17_GL000205v2_random` `chr17_GL000258v2_alt` `chr17_GL383563v3_alt` `chr17_GL383564v2_alt` `chr17_JH159146v1_alt` `chr17_JH159147v1_alt` `chr17_JH159148v1_alt` `chr17_KI270857v1_alt` `chr17_KI270858v1_alt` `chr17_KI270860v1_alt` `chr17_KI270861v1_alt` `chr17_KI270862v1_alt` `chr17_KI270907v1_alt` `chr17_KI270908v1_alt` `chr17_KI270909v1_alt` `chr17_KI270910v1_alt` `chr18` `chr18_GL383567v1_alt` `chr18_GL383569v1_alt` `chr18_GL383571v1_alt` `chr18_GL383572v1_alt` `chr18_KI270911v1_alt` `chr18_KI270912v1_alt` `chr19` `chr19_GL000209v2_alt` `chr19_GL383573v1_alt` `chr19_GL383574v1_alt` `chr19_GL383575v2_alt` `chr19_GL383576v1_alt` `chr19_GL949746v1_alt` `chr19_GL949747v2_alt` `chr19_GL949748v2_alt` `chr19_GL949749v2_alt` `chr19_GL949750v2_alt` `chr19_GL949751v2_alt` `chr19_GL949752v1_alt` `chr19_GL949753v2_alt` `chr19_KI270865v1_alt` `chr19_KI270867v1_alt` `chr19_KI270868v1_alt` `chr19_KI270882v1_alt` `chr19_KI270883v1_alt` `chr19_KI270884v1_alt` `chr19_KI270885v1_alt` `chr19_KI270886v1_alt` `chr19_KI270887v1_alt` `chr19_KI270888v1_alt` `chr19_KI270889v1_alt` `chr19_KI270890v1_alt` `chr19_KI270891v1_alt` `chr19_KI270914v1_alt` `chr19_KI270915v1_alt` `chr19_KI270916v1_alt` `chr19_KI270917v1_alt` `chr19_KI270918v1_alt` `chr19_KI270919v1_alt` `chr19_KI270920v1_alt` `chr19_KI270921v1_alt` `chr19_KI270922v1_alt` `chr19_KI270923v1_alt` `chr19_KI270929v1_alt` `chr19_KI270930v1_alt` `chr19_KI270931v1_alt` `chr19_KI270932v1_alt` `chr19_KI270933v1_alt` `chr19_KI270938v1_alt` `chr1_GL383518v1_alt` `chr1_GL383519v1_alt` `chr1_GL383520v2_alt` `chr1_KI270706v1_random` `chr1_KI270707v1_random` `chr1_KI270708v1_random` `chr1_KI270711v1_random` `chr1_KI270713v1_random` `chr1_KI270714v1_random` `chr1_KI270759v1_alt` `chr1_KI270761v1_alt` `chr1_KI270762v1_alt` `chr1_KI270763v1_alt` `chr1_KI270765v1_alt` `chr1_KI270766v1_alt` `chr1_KI270892v1_alt` `chr2` `chr20` `chr20_GL383577v2_alt` `chr20_KI270869v1_alt` `chr20_KI270870v1_alt` `chr21` `chr21_GL383579v2_alt` `chr21_GL383580v2_alt` `chr21_GL383581v2_alt` `chr21_KI270872v1_alt` `chr22` `chr22_GL383582v2_alt` `chr22_GL383583v2_alt` `chr22_KB663609v1_alt` `chr22_KI270731v1_random` `chr22_KI270733v1_random` `chr22_KI270734v1_random` `chr22_KI270875v1_alt` `chr22_KI270876v1_alt` `chr22_KI270877v1_alt` `chr22_KI270878v1_alt` `chr22_KI270879v1_alt` `chr22_KI270928v1_alt` `chr2_GL383521v1_alt` `chr2_GL383522v1_alt` `chr2_GL582966v2_alt` `chr2_KI270768v1_alt` `chr2_KI270769v1_alt` `chr2_KI270774v1_alt` `chr2_KI270776v1_alt` `chr2_KI270894v1_alt` `chr3` `chr3_GL383526v1_alt` `chr3_JH636055v2_alt` `chr3_KI270777v1_alt` `chr3_KI270779v1_alt` `chr3_KI270780v1_alt` `chr3_KI270782v1_alt` `chr3_KI270784v1_alt` `chr3_KI270895v1_alt` `chr3_KI270924v1_alt` `chr3_KI270934v1_alt` `chr3_KI270935v1_alt` `chr3_KI270936v1_alt` `chr3_KI270937v1_alt` `chr4` `chr4_GL000008v2_random` `chr4_GL000257v2_alt` `chr4_GL383527v1_alt` `chr4_GL383528v1_alt` `chr4_KI270786v1_alt` `chr4_KI270789v1_alt` `chr4_KI270896v1_alt` `chr4_KI270925v1_alt` `chr5` `chr5_GL339449v2_alt` `chr5_GL383531v1_alt` `chr5_GL383532v1_alt` `chr5_GL949742v1_alt` `chr5_KI270791v1_alt` `chr5_KI270792v1_alt` `chr5_KI270793v1_alt` `chr5_KI270794v1_alt` `chr5_KI270795v1_alt` `chr5_KI270796v1_alt` `chr5_KI270897v1_alt` `chr5_KI270898v1_alt` `chr6` `chr6_GL000250v2_alt` `chr6_GL000251v2_alt` `chr6_GL000252v2_alt` `chr6_GL000253v2_alt` `chr6_GL000254v2_alt` `chr6_GL000255v2_alt` `chr6_GL000256v2_alt` `chr6_GL383533v1_alt` `chr6_KI270758v1_alt` `chr6_KI270797v1_alt` `chr6_KI270798v1_alt` `chr6_KI270800v1_alt` `chr6_KI270801v1_alt` `chr7` `chr7_GL383534v2_alt` `chr7_KI270803v1_alt` `chr7_KI270806v1_alt` `chr7_KI270808v1_alt` `chr7_KI270809v1_alt` `chr8` `chr8_KI270810v1_alt` `chr8_KI270811v1_alt` `chr8_KI270813v1_alt` `chr8_KI270814v1_alt` `chr8_KI270815v1_alt` `chr8_KI270816v1_alt` `chr8_KI270817v1_alt` `chr8_KI270819v1_alt` `chr8_KI270821v1_alt` `chr8_KI270822v1_alt` `chr8_KI270900v1_alt` `chr8_KI270926v1_alt` `chr9` `chr9_GL383539v1_alt` `chr9_GL383540v1_alt` `chr9_GL383541v1_alt` `chr9_GL383542v1_alt` `chr9_KI270823v1_alt` `chrM` `chrUn_GL000195v1` `chrUn_GL000213v1` `chrUn_GL000216v2` `chrUn_GL000218v1` `chrUn_GL000219v1` `chrUn_GL000220v1` `chrUn_GL000224v1` `chrUn_KI270442v1` `chrUn_KI270741v1` `chrUn_KI270743v1` `chrUn_KI270744v1` `chrUn_KI270750v1` `chrUn_KI270752v1` `chrX` `chrX_KI270880v1_alt` `chrX_KI270881v1_alt` `chrX_KI270913v1_alt` `chrY`

### Prefix annotations

NA

### A transcript on positive strand

```text
chr1	ensGene	transcript	173477266	173488807	.	+	.	gene_id "ENSG00000117592"; transcript_id "ENST00000340385";  gene_name "ENSG00000117592";
chr1	ensGene	exon	173477266	173477492	.	+	.	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "1"; exon_id "ENST00000340385.1"; gene_name "ENSG00000117592";
chr1	ensGene	5UTR	173477266	173477397	.	+	.	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "1"; exon_id "ENST00000340385.1"; gene_name "ENSG00000117592";
chr1	ensGene	CDS	173477398	173477492	.	+	0	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "1"; exon_id "ENST00000340385.1"; gene_name "ENSG00000117592";
chr1	ensGene	exon	173481326	173481482	.	+	.	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "2"; exon_id "ENST00000340385.2"; gene_name "ENSG00000117592";
chr1	ensGene	CDS	173481326	173481482	.	+	1	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "2"; exon_id "ENST00000340385.2"; gene_name "ENSG00000117592";
chr1	ensGene	exon	173485361	173485507	.	+	.	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "3"; exon_id "ENST00000340385.3"; gene_name "ENSG00000117592";
chr1	ensGene	CDS	173485361	173485507	.	+	0	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "3"; exon_id "ENST00000340385.3"; gene_name "ENSG00000117592";
chr1	ensGene	exon	173486255	173486401	.	+	.	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "4"; exon_id "ENST00000340385.4"; gene_name "ENSG00000117592";
chr1	ensGene	CDS	173486255	173486401	.	+	0	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "4"; exon_id "ENST00000340385.4"; gene_name "ENSG00000117592";
chr1	ensGene	exon	173487735	173488807	.	+	.	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "5"; exon_id "ENST00000340385.5"; gene_name "ENSG00000117592";
chr1	ensGene	CDS	173487735	173487860	.	+	0	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "5"; exon_id "ENST00000340385.5"; gene_name "ENSG00000117592";
chr1	ensGene	3UTR	173487864	173488807	.	+	.	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "5"; exon_id "ENST00000340385.5"; gene_name "ENSG00000117592";
chr1	ensGene	start_codon	173477398	173477400	.	+	0	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "1"; exon_id "ENST00000340385.1"; gene_name "ENSG00000117592";
chr1	ensGene	stop_codon	173487861	173487863	.	+	0	gene_id "ENSG00000117592"; transcript_id "ENST00000340385"; exon_number "5"; exon_id "ENST00000340385.5"; gene_name "ENSG00000117592";
```

### A transcript on negative strand

```text
chr1	ensGene	transcript	14404	29570	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147";  gene_name "ENSG00000227232";
chr1	ensGene	exon	14404	14501	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "1"; exon_id "ENST00000488147.1"; gene_name "ENSG00000227232";
chr1	ensGene	exon	15005	15038	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "2"; exon_id "ENST00000488147.2"; gene_name "ENSG00000227232";
chr1	ensGene	exon	15796	15947	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "3"; exon_id "ENST00000488147.3"; gene_name "ENSG00000227232";
chr1	ensGene	exon	16607	16765	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "4"; exon_id "ENST00000488147.4"; gene_name "ENSG00000227232";
chr1	ensGene	exon	16858	17055	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "5"; exon_id "ENST00000488147.5"; gene_name "ENSG00000227232";
chr1	ensGene	exon	17233	17368	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "6"; exon_id "ENST00000488147.6"; gene_name "ENSG00000227232";
chr1	ensGene	exon	17606	17742	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "7"; exon_id "ENST00000488147.7"; gene_name "ENSG00000227232";
chr1	ensGene	exon	17915	18061	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "8"; exon_id "ENST00000488147.8"; gene_name "ENSG00000227232";
chr1	ensGene	exon	18268	18366	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "9"; exon_id "ENST00000488147.9"; gene_name "ENSG00000227232";
chr1	ensGene	exon	24738	24891	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "10"; exon_id "ENST00000488147.10"; gene_name "ENSG00000227232";
chr1	ensGene	exon	29534	29570	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000488147"; exon_number "11"; exon_id "ENST00000488147.11"; gene_name "ENSG00000227232";
```
