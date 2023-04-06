all: \
	GCF_000001405.40_GRCh38.p14_genomic.gtf \
	GCA_000001405.29_GRCh38.p14_genomic.gtf \
	gencode.v43.primary_assembly.annotation.gtf \
	gencode.v43.chr_patch_hapl_scaff.annotation.gtf \
	gencode.v43.primary_assembly.basic.annotation.gtf  \
	gencode.v43.chr_patch_hapl_scaff.basic.annotation.gtf \
	hg38.ensGene.gtf \
	hg38.knownGene.gtf \
	hg38.ncbiRefSeq.gtf \
	hg38.refGene.gtf \
	Homo_sapiens.GRCh38.109.gtf \
	GRCh38_GENCODE_rmsk_TE.gtf \
	GRCh38_Ensembl_rmsk_TE.gtf \
	hg38.fa \
	hg38.masked.fa \
	GRCh38_full_analysis_set_plus_decoy_hla.fa \
	Homo_sapiens.GRCh38.dna.toplevel.fa \
	Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	Homo_sapiens.GRCh38.dna_sm.toplevel.fa \
	Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
	Homo_sapiens.GRCh38.dna_rm.toplevel.fa \
	Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa \
	GRCh38.p13.genome.fa \
	GCF_000001405.40_GRCh38.p14_genomic.fna \
	GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna \
	GCA_000001405.29_GRCh38.p14_genomic.fna \
	hg38.analysisSet.fa

%: %.gz
	gunzip $<

GCF_000001405.40_GRCh38.p14_genomic.gtf.gz:
	axel "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"

GCA_000001405.29_GRCh38.p14_genomic.gtf.gz:
	axel "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.gtf.gz"

gencode.v43.primary_assembly.annotation.gtf.gz:
	axel "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz"

gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz:
	axel "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz"

gencode.v43.primary_assembly.basic.annotation.gtf.gz:
	axel "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.basic.annotation.gtf.gz"

gencode.v43.chr_patch_hapl_scaff.basic.annotation.gtf.gz:
	axel "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.chr_patch_hapl_scaff.basic.annotation.gtf.gz"


hg38.ensGene.gtf.gz:
	axel "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz"

hg38.knownGene.gtf.gz:
	axel "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz"

hg38.ncbiRefSeq.gtf.gz:
	axel "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz"

hg38.refGene.gtf.gz:
	axel "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz"

Homo_sapiens.GRCh38.109.gtf.gz:
	axel "http://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"

GRCh38_GENCODE_rmsk_TE.gtf.gz:
	axel "http://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCh38_GENCODE_rmsk_TE.gtf.gz"

GRCh38_Ensembl_rmsk_TE.gtf.gz:
	axel "https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCh38_Ensembl_rmsk_TE.gtf.gz"

hg38.fa.gz:
	axel "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz"

hg38.fa.masked.gz:
	axel "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.masked.gz"

hg38.masked.fa: hg38.fa.masked
	mv hg38.fa.masked hg38.masked.fa

Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz:
	axel "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz"

Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz:
	axel "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"

Homo_sapiens.GRCh38.dna_rm.toplevel.fa.gz:
	axel "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.fa.gz"

Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.gz:
	axel "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.gz"

Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz:
	axel "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

Homo_sapiens.GRCh38.dna.toplevel.fa.gz:
	axel "ttps://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"

GRCh38.p13.genome.fa.gz:
	axel "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.p13.genome.fa.gz"

GRCh38.primary_assembly.genome.fa.gz:
	axel "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz"

GCF_000001405.40_GRCh38.p14_genomic.fna.gz:
	axel "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"

GCA_000001405.29_GRCh38.p14_genomic.fna.gz:
	axel "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz"

GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz:
	axel "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz"

hg38.analysisSet.fa.gz:
	axel "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz"

GRCh38_full_analysis_set_plus_decoy_hla.fa:
	axel "https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
