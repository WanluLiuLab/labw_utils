#!/usr/bin/env bash
rm -f bench_result.tsv
python bench.py
Rscript plot.R
