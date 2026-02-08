#!/bin/bash -e

cd ../big_hits

python ../py/make_hist_fold_sf_and_scop40c.py \
	--hits ../big_data/tm.scop40 \
	--type tm \
	--fields 1 2 3 \
	--output ../hist/tm.scop40

python ../py/make_hist_fold_sf_and_scop40c.py \
	--hits foldseek.scop40 \
	--type bits \
	--fields 1 2 3 \
	--output ../hist/foldseek_bits.scop40

python ../py/make_hist_fold_sf_and_scop40c.py \
	--hits blastp.scop40 \
	--type blastp \
	--fields 1 2 4 \
	--output ../hist/blast_rawscore.scop40

python ../py/make_hist_fold_sf_and_scop40c.py \
	--hits blastp.scop40 \
	--type blastp \
	--fields 1 2 4 \
	--minlength 150 \
	--maxlength 250 \
	--output ../hist/blast_rawscore.minlength150.maxlength250.scop40

python ../py/make_hist_fold_sf_and_scop40c.py \
    --hits reseek27_verysensitive.scop40 \
    --type ts \
    --fields 1 2 3 \
    --output ../hist/reseek_ts.scop40

python ../py/make_hist_fold_sf_and_scop40c.py \
    --hits ../data/dali.scop40 \
    --type dali \
    --fields 1 2 3 \
    --output ../hist/dali.scop40
