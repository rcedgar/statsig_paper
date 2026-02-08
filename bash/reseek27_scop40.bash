#!/bin/bash -e

mkdir -p ../big_hits
cd ../big_hits

$src/reseek/github_releases/reseek-v2.7-linux-x86 \
	-search ../data/scop40.bca \
	-db ../data/scop40.bca \
	-verysensitive \
	-columns query+target+newts \
	-output reseek27_verysensitive.scop40

$src/reseek/github_releases/reseek-v2.7-linux-x86 \
	-search ../data/scop40.bca \
	-db ../data/scop40.bca \
	-fast \
	-columns query+target+evalue+pvalue \
	-output reseek27_fast.scop40
