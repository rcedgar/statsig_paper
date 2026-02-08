#!/bin/bash -e

mkdir -p ../big_hits

blastp \
	-query ../data/scop40.fa \
	-subject ../data/scop40.fa \
	-evalue 10 \
	-outfmt "6 qseqid sseqid evalue score" \
	> ../big_hits/blastp.scop40

ls -lh ../big_hits/blastp.scop40
