#!/bin/bash -e

foldseek=foldseek-v10
tmpdir=../tmp_foldseek
mkdir -p $tmpdir

mkdir -p ../big_hits

$foldseek \
	  easy-search \
	  ../big_foldseek/scop40/scop40db \
	  ../big_foldseek/scop40/scop40db \
	  ../big_hits/foldseek.scop40 \
	  $tmpdir \
	  -e 1e9 \
	  --format-output query,target,bits,evalue \
	  --exhaustive-search

rm -rf $tmpdir

$foldseek \
	  easy-search \
	  ../big_foldseek/scop40/scop40db \
	  ../big_foldseek/scop40/scop40db \
	  ../big_hits/foldseek.scop40.prob \
	  $tmpdir \
	  --format-output query,target,prob,evalue,bits

rm -rf $tmpdir
