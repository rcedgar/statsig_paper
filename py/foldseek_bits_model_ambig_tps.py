#!/usr/bin/python3
import sys
import numpy as np
from lookup import read_lookup
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')

hitsfn = "../big_hits/foldseek.scop40.bits_ge40"
lookupfn = "../data/scop40c.lookup"

plotfn = "../figures/foldseek_bits_model_ambig_tps.svg"
tsvfn = "../results/foldseek_bits_model_ambig_tps.tsv"

ftsv = open(tsvfn, "w")

plotw = 4
ploth = 4

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(plotw, ploth))

dom2sf, dom2fold = read_lookup(lookupfn)

maxbitscore = 100
minbitscore = 40

color_with_discards = "black"
color_no_discards = "black"

lw = 2
ls_with_discards = "dotted"
ls_no_discards = "solid"
ls_model = "dashed"

model_names =	[ "Half ambig.\nare FP",	"Weighted\nmodel"	]
model_colors =	[ "gray",	"gray"	]
model_ls =		[ "dotted",	"dashed"]
p1s =			[ 0.5,		0.6		]
p2s =			[ 0.5,		0.2		]

nmodels = len(model_names)

def get_prob_diff_same_fold_is_fp_model(bitscore, model_idx):
	p1 = p1s[model_idx]
	p2 = p2s[model_idx]
	r = (bitscore - minbitscore)/(maxbitscore - minbitscore)
	assert r >= 0 and r <= 1
	prob = (1 - r)*p1 + r*p2
	assert prob <= p1 and prob >= p2
	return prob

counts_same_sf = [0]*(maxbitscore+1)
counts_diff_fold = [0]*(maxbitscore+1)
counts_diff_sf_same_fold = [0]*(maxbitscore+1)

for line in open(hitsfn):
	q, t, sbitscore = line[:-1].split('\t')
	if q == t:
		continue
	qsf = dom2sf.get(q)
	tsf = dom2sf.get(t)
	if qsf is None or tsf is None:
		continue

	qfold = dom2fold.get(q)
	tfold = dom2fold.get(t)
	bitscore = int(sbitscore)
	if bitscore > maxbitscore:
		continue

	if qsf == tsf:
		counts_same_sf[bitscore] += 1
	elif qfold != tfold:
		counts_diff_fold[bitscore] += 1
	elif qfold == tfold:
		counts_diff_sf_same_fold[bitscore] += 1
	else:
		assert False

s = "bits"
s += "\tprob_fp_no_discards"
s += "\tprob_fp_after_discards"

for model_idx in range(nmodels):
	s += "\tP%d" % model_idx

for model_idx in range(nmodels):
	s += "\tprob_fp%d" % model_idx
ftsv.write(s + '\n')

bitscores = []

probs_fp_no_discards = []
probs_fp_after_discards = []
probs_fp_model_vec = []
for model_idx in range(nmodels):
	probs_fp_model_vec.append([])

for bitscore in range(maxbitscore, minbitscore-1, -1):
	bitscores.append(bitscore)

	n_same_sf = counts_same_sf[bitscore]
	n_diff_fold = counts_diff_fold[bitscore]
	n_diff_sf_same_fold = counts_diff_sf_same_fold[bitscore]

	n_after_discards = n_same_sf + n_diff_fold
	prob_fp_after_discards = n_diff_fold/n_after_discards

	n_diff_sf = n_diff_fold + n_diff_sf_same_fold
	n = n_same_sf + n_diff_fold + n_diff_sf_same_fold
	prob_fp_no_discards = n_diff_sf/n

	probs_fp_no_discards.append(prob_fp_no_discards)
	probs_fp_after_discards.append(prob_fp_after_discards)

	s = "%d" % bitscore
	s += "\t%.3g" % prob_fp_no_discards
	s += "\t%.3g" % prob_fp_after_discards

	for model_idx in range(nmodels):
		model_prob = get_prob_diff_same_fold_is_fp_model(bitscore, model_idx)
		s += "\t%.4f" % model_prob

	for model_idx in range(nmodels):
		model_prob = get_prob_diff_same_fold_is_fp_model(bitscore, model_idx)
		n_diff_sf_model = n_diff_fold + n_diff_sf_same_fold*model_prob
		n_diff_fold_model = n_diff_fold + n_diff_sf_same_fold*(1 - model_prob)
		prob_fp_model = n_diff_sf_model/n
		probs_fp_model_vec[model_idx].append(prob_fp_model)

		s += "\t%.3g" % prob_fp_model
	ftsv.write(s + '\n')

print("Writing " + tsvfn)
ftsv.close()

ax.plot(bitscores, probs_fp_no_discards, label="Keep ambig.", 
		lw=lw, color=color_no_discards, ls=ls_no_discards)

for model_idx in range(nmodels):
	ax.plot(bitscores, probs_fp_model_vec[model_idx], lw=lw, ls=model_ls[model_idx],
		color=model_colors[model_idx], label=model_names[model_idx])

ax.plot(bitscores, probs_fp_after_discards, label="Discard ambig.\n(Foldseek)", 
		lw=lw, color=color_with_discards, ls=ls_with_discards)

ax.legend()
ax.set_xlabel("Bits")
ax.set_ylabel("P(FP)")

plt.tight_layout()
print("Writing " + plotfn)
plt.savefig(plotfn)
