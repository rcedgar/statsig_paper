#!/bin/bash -e

mkdir -p ../big_hits

reseek=../bin/reseek_7ba344c

$reseek \
	-para_scop40 ../data/scop40.3di.fa \
	-mxname 3Di \
	-alignmethod para \
	-output ../big_hits/3di.scop40 \
	-seqsmethod 3Di

$reseek \
	-para_scop40 ../data/scop40.3di.shuffled.fa \
	-mxname 3Di \
	-alignmethod para \
	-output ../big_hits/3di.scop40.shuffled \
	-seqsmethod 3Di

$reseek \
	-para_scop40 ../data/scop40.3di.shuffled_all.fa \
	-mxname 3Di \
	-alignmethod para \
	-output ../big_hits/3di.scop40.shuffled_all \
	-seqsmethod 3Di

$reseek \
	-para_scop40 ../data/scop40.3di.shuffled_lc.fa \
	-mxname 3Di \
	-alignmethod para \
	-output ../big_hits/3di.scop40.shuffled_lc \
	-seqsmethod 3Di

$reseek \
	-para_scop40 ../data/scop40.3di.shuffled1k.fa \
	-mxname 3Di \
	-alignmethod para \
	-output ../big_hits/3di.scop40.shuffled1k \
	-seqsmethod 3Di
