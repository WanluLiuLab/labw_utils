#!/usr/bin/env bash

if [ ! -f L4_rep2.bam ]; then
    {
        if [ ! -f ce11.ncbiRefSeq.chr1.gtf ]; then
            axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz &>> /dev/null
            gunzip ce11.ncbiRefSeq.gtf.gz
        fi
        if [ ! -f ce11.chr1.fa ]; then
            axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz &>> /dev/null
            gunzip ce11.fa.gz
        fi
        axel ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3245471/L4_rep2.fastq.gz &>> /dev/null
        gunzip L4_rep2.fastq.gz
        minimap2 -a -x splice ce11.fa L4_rep2.fastq | \
        samtools sort -o L4_rep2.bam
        samtools index L4_rep2.bam
     }  &>> preparation.log
fi

if [ ! -d L4_rep2.fastq.stats.d ]; then
    python -m labw_utils.bioutils describe_fastq L4_rep2.fastq &>> describe_fastq.log
fi
if [ ! -f get_exonic_depth.log ]; then 
    python -m labw_utils.bioutils get_exonic_depth \
        -s L4_rep2.bam \
        -g ce11.ncbiRefSeq.gtf | {
    grep -v 'inferred from feature transcript' || true
} &>> get_exonic_depth.log
fi
