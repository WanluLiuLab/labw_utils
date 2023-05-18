import os

from pyspark.sql import SparkSession

from labw_utils.commonutils.libfrontend import setup_basic_logger
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger()

PRIMIARY_KEY_RENAME_TABLE = {
    "transcript_id": "ensdb_pk_transcript_id",
    "gene_id": "ensdb_pk_gene_id",
    "exon_id": "ensdb_pk_exon_id",
    "repeat_feature_id": "ensdb_pk_repeat_feature_id"
}

setup_basic_logger()
spark = SparkSession.builder.getOrCreate()
_lh.info("Using pyspark version %s", spark.version)
transcript_table = spark.read.parquet(os.path.join("converted_parquet", "transcript.parquet"))
seq_region_synonym_table = spark.read.parquet(os.path.join("converted_parquet", "seq_region_synonym.parquet"))
external_db_table = spark.read.parquet(os.path.join("converted_parquet", "external_db.parquet"))
gene_table = spark.read.parquet(os.path.join("converted_parquet", "gene.parquet"))
exon_table = spark.read.parquet(os.path.join("converted_parquet", "exon.parquet"))
exon_transcript_table = spark.read.parquet(os.path.join("converted_parquet", "exon_transcript.parquet"))
repeat_feature_table = spark.read.parquet(os.path.join("converted_parquet", "repeat_feature.parquet"))
repeat_consensus_table = spark.read.parquet(os.path.join("converted_parquet", "repeat_consensus.parquet"))

seq_region_synonym_converter = (
    seq_region_synonym_table.
    select(
        "seq_region_id",
        "external_db_id",
        "synonym"
    ).
    join(
        external_db_table.
        select("external_db_id", "db_name").
        filter("db_name = 'UCSC'"),
        on="external_db_id"
    ).
    select(
        "seq_region_id",
        "synonym"
    )
)


def merge_gene():
    global ens_hgnc_gene_map_table
    _lh.info("Merging genes...")
    final_locus = (
        gene_table.
        select(
            "gene_id",
            "seq_region_id",
            "seq_region_start",
            "seq_region_end",
            "seq_region_strand"
        ).
        join(
            seq_region_synonym_converter,
            on="seq_region_id",
            how="inner"
        ).
        drop("seq_region_id").
        withColumnRenamed("synonym", "seqname").
        withColumnRenamed("seq_region_start", "start").
        withColumnRenamed("seq_region_end", "end").
        withColumnRenamed("seq_region_strand", "strand").
        withColumnsRenamed(PRIMIARY_KEY_RENAME_TABLE)
    )
    final_gene = (
        gene_table.
        select(
            "gene_id",
            "stable_id",
            "version",
            "created_date",
            "modified_date"
        ).
        withColumnRenamed("stable_id", "ensdb_gene_id").
        withColumnRenamed("version", "ensdb_gene_version").
        withColumnRenamed("created_date", "ensdb_gene_created_date").
        withColumnRenamed("modified_date", "ensdb_gene_modified_date").
        withColumnRenamed("gene_id", "ensdb_pk_gene_id")
    )
    joint_table = final_gene.join(final_locus, on="ensdb_pk_gene_id", how="inner")
    joint_table.write.parquet(
        "ensdb_genes.parquet.d",
        mode="overwrite"
    )
    return joint_table


def merge_transcripts():
    _lh.info("Merging transcripts...")
    final_locus = (
        transcript_table.
        select(
            "transcript_id",
            "seq_region_id",
            "seq_region_start",
            "seq_region_end",
            "seq_region_strand"
        ).
        join(
            seq_region_synonym_converter,
            on="seq_region_id",
            how="inner"
        ).
        drop("seq_region_id").
        withColumnRenamed("synonym", "seqname").
        withColumnRenamed("seq_region_start", "start").
        withColumnRenamed("seq_region_end", "end").
        withColumnRenamed("seq_region_strand", "strand").
        withColumnsRenamed(PRIMIARY_KEY_RENAME_TABLE)
    )
    final_misc = (
        transcript_table.
        select(
            "transcript_id",
            "gene_id",
            "source",
            "description",
            "stable_id",
            "version",
            "created_date",
            "modified_date"
        ).
        withColumnRenamed("stable_id", "ensdb_transcript_id").
        withColumnRenamed("version", "ensdb_transcript_version").
        withColumnRenamed("created_date", "ensdb_transcript_created_date").
        withColumnRenamed("modified_date", "ensdb_transcript_modified_date").
        withColumnsRenamed(PRIMIARY_KEY_RENAME_TABLE)
    )
    joint_table = final_locus.join(final_misc, on="ensdb_pk_transcript_id", how="inner")
    joint_table.write.parquet(
        "ensdb_transcripts.parquet.d",
        mode="overwrite"
    )
    return joint_table


def merge_exons():
    _lh.info("Merging exons...")
    final_locus = (
        exon_table.
        select(
            "exon_id",
            "seq_region_id",
            "seq_region_start",
            "seq_region_end",
            "seq_region_strand"
        ).
        join(
            seq_region_synonym_converter,
            on="seq_region_id",
            how="inner"
        ).
        drop("seq_region_id").
        withColumnRenamed("synonym", "seqname").
        withColumnRenamed("seq_region_start", "start").
        withColumnRenamed("seq_region_end", "end").
        withColumnRenamed("seq_region_strand", "strand").
        withColumnsRenamed(PRIMIARY_KEY_RENAME_TABLE)
    )
    final_misc = (
        exon_table.
        select(
            "exon_id",
            "stable_id",
            "version",
            "created_date",
            "modified_date"
        ).
        withColumnRenamed("stable_id", "ensdb_exon_id").
        withColumnRenamed("version", "ensdb_exon_version").
        withColumnRenamed("created_date", "ensdb_exon_created_date").
        withColumnRenamed("modified_date", "ensdb_exon_modified_date").
        withColumnsRenamed(PRIMIARY_KEY_RENAME_TABLE)
    )
    final_mapping = (
        exon_transcript_table.
        select(
            "transcript_id",
            "exon_id"
        ).
        withColumnsRenamed(PRIMIARY_KEY_RENAME_TABLE)
    )

    joint_table = (
        final_locus.
        join(final_mapping, on="ensdb_pk_exon_id", how="inner").
        join(final_misc, on="ensdb_pk_exon_id", how="inner")
    )
    joint_table.write.parquet(
            "ensdb_exons.parquet.d",
            mode="overwrite"
        )
    return joint_table


def merge_repeats():
    _lh.info("Merging repeats...")
    final_locus = (
        repeat_feature_table.
        select(
            "repeat_feature_id",
            "seq_region_id",
            "seq_region_start",
            "seq_region_end",
            "seq_region_strand"
        ).
        join(
            seq_region_synonym_converter,
            on="seq_region_id",
            how="inner"
        ).
        drop("seq_region_id").
        withColumnRenamed("synonym", "seqname").
        withColumnRenamed("seq_region_start", "start").
        withColumnRenamed("seq_region_end", "end").
        withColumnRenamed("seq_region_strand", "strand").
        withColumnsRenamed(PRIMIARY_KEY_RENAME_TABLE)
    )
    final_misc = (
        repeat_feature_table.
        select(
            "repeat_feature_id",
            "repeat_consensus_id",
            "score",
            "repeat_start",
            "repeat_end"
        ).
        join(
            repeat_consensus_table.drop("repeat_consensus"),
            on="repeat_consensus_id",
            how="inner"
        ).
        withColumnsRenamed(PRIMIARY_KEY_RENAME_TABLE)
    )
    joint_table = final_locus.join(final_misc, on="ensdb_pk_repeat_feature_id", how="inner")
    joint_table.write.parquet(
            "ensdb_repeats.parquet.d",
            mode="overwrite"
        )
    return joint_table



if __name__ == "__main__":
    gene_table = merge_gene()
    transcripts_table = merge_transcripts()
    exons_table = merge_exons()
    repeats_table = merge_repeats()

