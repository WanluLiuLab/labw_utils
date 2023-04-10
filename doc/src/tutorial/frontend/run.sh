#!/usr/bin/env bash
set -ue
if [ ! -f ce11.ncbiRefSeq.gtf ]; then
    axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz &>>/dev/null
    gunzip -f ce11.ncbiRefSeq.gtf.gz
fi
if [ ! -f ce11.fa ]; then
    axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz &>>/dev/null
    gunzip -f ce11.fa.gz
fi
if [ ! -f L4_rep2.bam ]; then
    {
        axel ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3245471/L4_rep2.fastq.gz &>>/dev/null
        gunzip -f L4_rep2.fastq.gz
        minimap2 -a -x splice ce11.fa L4_rep2.fastq |
            samtools sort -o L4_rep2.bam
        samtools index L4_rep2.bam
    } &>>preparation.log
fi

if [ ! -d L4_rep2.fastq.stats.d ]; then
    python -m labw_utils.bioutils describe_fastq L4_rep2.fastq &>>describe_fastq.log
fi

if [ ! -f get_exonic_depth.finished ]; then
    python -m labw_utils.bioutils get_exonic_depth \
        -s L4_rep2.bam \
        -g ce11.ncbiRefSeq.gtf 2>&1 | {
        grep -v 'inferred from feature transcript' || true
    } &>>get_exonic_depth.log
    touch get_exonic_depth.finished
fi
if [ ! -f filter_gtf_by_attribute.finished ]; then
    python -m labw_utils.bioutils filter_gtf_by_attribute \
    -g ce11.ncbiRefSeq.gtf \
    --attribute_name gene_id \
    --attribute_values filter.regex \
    --out ce11.filtered.gtf \
    --regex &>>filter_gtf_by_attribute.log
    touch filter_gtf_by_attribute.finished
fi
if [ ! -f ce11_trans_filtered.fa ]; then
   python -m labw_utils.bioutils transcribe \
    -f ce11.fa \
    -g ce11.filtered.gtf \
    --no_write_single_transcript \
    -o ce11_trans_filtered.fa &>>transcribe.log
    touch ce11_trans_filtered.fa
fi
