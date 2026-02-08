#!/bin/bash -e

python ../py/hist_3di_selected_doms.py

for scop_class in a b c d ; do
	python ../py/hist_3di_selected_doms.py $scop_class
done
