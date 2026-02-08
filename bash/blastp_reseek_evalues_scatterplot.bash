#!/bin/bash -e

python ../py/blastp_reseek_evalues_sample.py \
	> ../results//blastp_reseek_evalues_sample.tsv

python ../py/blastp_reseek_evalue_scatterplot.py
