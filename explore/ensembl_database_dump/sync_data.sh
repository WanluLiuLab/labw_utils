#!/usr/bin/env bash

rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_mysql/homo_sapiens_core_109_38 rsync_data
for fn in rsync_data/*.gz; do gunzip "${fn}" & done
wait

Rscript get_gene_ens_hgnc_mapping.R

docker run \
    --detach \
    -p 13306:3306 \
    --name labw_utils_ensdb \
    --env MYSQL_USER=ensdb \
    --env MYSQL_PASSWORD=ensdb \
    --env MYSQL_ROOT_PASSWORD=ensdb \
    --env MYSQL_DATABASE=homo_sapiens_core_109_38 \
    mariadb:10.0.30

mariadb \
    --host=localhost \
    --port=13306 \
    --password=ensdb \
    --database=homo_sapiens_core_109_38 \
    --user=ensdb \
    < rsync_data/homo_sapiens_core_109_38.sql
mysqlimport \
    --host=localhost \
    --port=13306 \
    --password=ensdb \
    --user=ensdb \
    --fields-terminated-by='\t' \
    --fields-escaped-by='\\' \
    homo_sapiens_core_109_38 \
    -L rsync_data/*.txt
python convert_data_from_mysql_to_parquet.py

# shellcheck disable=SC2155
export SPARK_CONF_DIR="$(pwd)"
spark-submit \
    --executor-memory 5G \
    --driver-memory 20G \
    --files "$(pwd)/log4j.properties" \
    --conf spark.driver.extraJavaOptions="-Dlog4j.configuration=file://$(pwd)/log4j.properties" \
    --conf spark.executor.extraJavaOptions="-Dlog4j.configuration=file://$(pwd)/log4j.properties" \
    --conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
    --packages com.esotericsoftware:kryo:4.0.2 \
    src/fitk/flitkdb/__init__.py


docker rm -f labw_utils_ensdb
