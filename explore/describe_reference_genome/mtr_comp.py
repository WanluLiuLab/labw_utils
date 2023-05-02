import functools
import gc
import glob
import itertools
import multiprocessing
import os
import json
import sys

import pandas as pd
from pyspark.sql import SparkSession, DataFrame
from pyspark.sql import types as st

from tqdm import tqdm
from labw_utils.commonutils.libfrontend import setup_basic_logger
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.shutil_helper import rm_rf
from labw_utils.typing_importer import Any, Iterable, Tuple

_lh = get_logger(__name__)
setup_basic_logger()
spark = SparkSession.builder.getOrCreate()
sc = spark.sparkContext
sc.setLogLevel("FATAL")
spark.conf.set("spark.sql.parquet.columnarReaderBatchSize", multiprocessing.cpu_count())

conf = {
    "SparkConf": dict(sc.getConf().getAll()),
    "ALL_ENV": dict(os.environ),
    "ENV": sc.environment,
    "SparkVersion": spark.version,
    "PyVersion": sys.version
}
with open("spark_conf.json", "w") as writer:
    json.dump(conf, writer, indent=4)

_lh.info("Using pyspark version %s", spark.version)


@functools.lru_cache(maxsize=0)
def read_parquet(fn: str) -> DataFrame:
    df = spark.read.parquet(*glob.glob(os.path.join(fn, "*.parquet")))
    cur_n_partitions = df.rdd.getNumPartitions()
    df = df.repartition(max(cur_n_partitions, 2 * multiprocessing.cpu_count()))
    return df


@functools.lru_cache(maxsize=0)
def cnt_parquet(fn: str) -> int:
    df = read_parquet(fn)
    fc = df.select("START").cache()
    reti = fc.count()
    del fc
    return reti


def compare_exon_def(
        fn1: str,
        fn2: str
):
    def select_rename(mtr: DataFrame, suff: int) -> DataFrame:
        sel_names = (
            "FTSHA",
            "EXON_BOUNDARIES",
            "START",
            "END",
            "STRAND",
            "SEQNAME",
            "SPLICE_SITES",
            "CDNA_BLAKE2B",
            "CDNA_UNSPLICED_BLAKE2B"
        )
        rename_dict = {
            sel_name: f"{sel_name}{suff}"
            for sel_name in sel_names
        }
        return (
            mtr.
            select(*sel_names).
            withColumnsRenamed(rename_dict)
        )

    _lh.info("Merging %s <-> %s", fn1, fn2)
    exon_def1 = select_rename(read_parquet(fn1), 1)
    exon_def2 = select_rename(read_parquet(fn2), 2)
    exon_def_merged = exon_def1.crossJoin(exon_def2)
    # print(exon_def_merged.schema)

    exon_def_merged = (
        exon_def_merged.
        withColumn(
            "REGIONAL_EQ",
            (
                    (exon_def_merged["START1"] == exon_def_merged["START2"]) &
                    (exon_def_merged["END1"] == exon_def_merged["END2"]) &
                    (exon_def_merged["STRAND1"] == exon_def_merged["STRAND2"]) &
                    (exon_def_merged["SEQNAME1"] == exon_def_merged["SEQNAME2"])
            )
        ).
        drop(
            "START1", "END1", "STRAND1", "SEQNAME1",
            "START2", "END2", "STRAND2", "SEQNAME2",
        )
    )
    exon_def_merged = (
        exon_def_merged.
        withColumn(
            "EXON_DEF_EQ",
            (
                    (exon_def_merged["EXON_BOUNDARIES1"] == exon_def_merged["EXON_BOUNDARIES2"]) &
                    exon_def_merged["REGIONAL_EQ"]
            )
        ).
        drop("EXON_BOUNDARIES1", "EXON_BOUNDARIES2").
        withColumn(
            "FTSHA_EQ",
            (
                    (exon_def_merged["FTSHA1"] == exon_def_merged["FTSHA2"])
            )
        ).
        drop("FTSHA1", "FTSHA2").
        drop("EXON_BOUNDARIES1", "EXON_BOUNDARIES2").
        withColumn(
            "SPLICE_SITE_EQ",
            (
                    (exon_def_merged["SPLICE_SITES1"] == exon_def_merged["SPLICE_SITES2"]) &
                    exon_def_merged["REGIONAL_EQ"]
            )
        ).
        drop("SPLICE_SITES1", "SPLICE_SITES2").
        withColumn(
            "CDNA_EQ",
            exon_def_merged["CDNA_BLAKE2B1"] == exon_def_merged["CDNA_BLAKE2B2"]
        ).
        drop("CDNA_BLAKE2B1", "CDNA_BLAKE2B2").
        withColumn(
            "CDNA_UNSPLICED_EQ",
            exon_def_merged["CDNA_UNSPLICED_BLAKE2B1"] == exon_def_merged["CDNA_UNSPLICED_BLAKE2B2"]
        ).
        drop("CDNA_UNSPLICED_BLAKE2B1", "CDNA_UNSPLICED_BLAKE2B2")
    )

    final_df = spark.createDataFrame(
        [(
            fn1,
            fn2,
            cnt_parquet(fn1),
            cnt_parquet(fn2),
            exon_def_merged.filter("FTSHA_EQ").select("FTSHA_EQ").cache().count(),
            exon_def_merged.filter("REGIONAL_EQ").select("REGIONAL_EQ").cache().count(),
            exon_def_merged.filter("EXON_DEF_EQ").select("EXON_DEF_EQ").cache().count(),
            exon_def_merged.filter("SPLICE_SITE_EQ").select("SPLICE_SITE_EQ").cache().count(),
            exon_def_merged.filter("CDNA_EQ").select("CDNA_EQ").cache().count(),
            exon_def_merged.filter("CDNA_UNSPLICED_EQ").select("CDNA_UNSPLICED_EQ").cache().count()
        )],
        schema=[
            "F1", "F2", "LF1", "LF2",
            "FTSHA_EQ",
            "REGIONAL_EQ",
            "EXON_DEF_EQ",
            "SPLICE_SITE_EQ",
            "CDNA_EQ",
            "CDNA_UNSPLICED_EQ"
        ]
    )
    del exon_def_merged
    gc.collect()
    return final_df


def product_with_no_overlap(it: Iterable[Any]) -> Iterable[Tuple[Any, Any]]:
    lit = list(it)
    for _p1 in range(0, len(lit)):
        for _p2 in range(_p1 + 1, len(lit)):
            yield lit[_p1], lit[_p2]


if __name__ == "__main__":
    files = []
    for _fn in glob.glob(os.path.join("pre_processed_mtr", "*.parquet.d")):
        files.append(_fn)
    _final_df = spark.createDataFrame(
        data=[],
        schema=st.StructType([
            st.StructField("F1", st.StringType(), False),
            st.StructField("F2", st.StringType(), False),
            st.StructField("LF1", st.IntegerType(), False),
            st.StructField("LF2", st.IntegerType(), False),
            st.StructField("FTSHA_EQ", st.IntegerType(), False),
            st.StructField("REGIONAL_EQ", st.IntegerType(), False),
            st.StructField("EXON_DEF_EQ", st.IntegerType(), False),
            st.StructField("SPLICE_SITE_EQ", st.IntegerType(), False),
            st.StructField("CDNA_EQ", st.IntegerType(), False),
            st.StructField("CDNA_UNSPLICED_EQ", st.IntegerType(), False)
        ]),
        verifySchema=False
    )
    rm_rf("summary_mtr_eq.parquet.d")
    for _fn1, _fn2 in tqdm(list(itertools.product(files, files))):
        compare_exon_def(_fn1, _fn2).write.parquet("summary_mtr_eq.parquet.d", mode="append")
    pd.read_parquet("summary_mtr_eq.parquet.d").to_parquet("summary_mtr_eq.parquet", index=False)
