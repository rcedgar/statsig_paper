import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats

doms = [ "d1itha_", "d1rrea_", "d3shua_" ]
titles = [ "d1itha_ (a.1) all-alpha", "d1rrea_ (d.14) alpha+beta", "d3shua_ (b.36) all-beta" ]
plotfn = "../figures/foldseek_cmdline_doms_scop_classes.svg"

ndoms = len(doms)
assert ndoms > 1

scop_class2fillcolor = { "a" : "lime", "b" : "orange", "all" : "lightgray" }
scop_class2color = { "a" : "lime", "b" : "orange", "all" : "gray" }
scop_class2name = { "a" : "all-alpha", "b" : "all-beta", "all" : "all classes" }
scop_class2ls = { "a" : "solid", "b" : "solid", "all" : "dashed" }

scop_classes = [ "a", "b", "all" ]

matplotlib.use('Agg')

plotw = 4*ndoms
ploth = 6

hitsfn = "../big_hits/foldseek.scop40"
lookupfn = "../data/scop40c.lookup"
scop_classes = [ "a", "b", "all" ]

dom2sf = {}
for line in open(lookupfn):
	xdom, scopid = line[:-1].split('\t')
	flds = scopid.split('.')
	sf = flds[0] + "." + flds[1] + "." + flds[2]
	dom2sf[xdom] = sf

minscore = -10
maxscore = 90
nbins = maxscore - minscore + 1
bin_mids = []
for score in range(minscore, maxscore+1):
	bin_mids.append(score)

dom2scop_class2counts = {}
for dom in doms:
	dom2scop_class2counts[dom] = {}
	for scop_class in scop_classes:
		dom2scop_class2counts[dom][scop_class] = [0]*nbins

for line in open(hitsfn):
	flds = line[:-1].split('\t')
	q = flds[0]
	if not q in doms:
		continue
	t = flds[1]
	qsf = dom2sf.get(q)
	tsf = dom2sf.get(t)
	if qsf is None or tsf is None or qsf == tsf:
		continue
	score = int(flds[2])
	if score < minscore or score > maxscore:
		continue
	binidx = score - minscore
	dom2scop_class2counts[q]["all"][binidx] += 1
	scop_class = tsf[0]
	if scop_class in scop_classes:
		dom2scop_class2counts[q][scop_class][binidx] += 1

def plot(ax, dom, scores, class2counts, xlo, xhi, ylo, yhi, title):
	for scop_class in scop_classes:
		counts = class2counts[scop_class]
		edgecolor = scop_class2color[scop_class]
		color = scop_class2color[scop_class]
		fill_color = scop_class2fillcolor[scop_class]
		name = scop_class2name[scop_class]
		ls = scop_class2ls[scop_class]

		ax.plot(scores, counts, color=color, ls=ls)
		if not fill_color is None:
			ax.fill_between(scores, counts, color=fill_color, alpha=0.2, ls=ls, label=name)
		ax.set_xlabel("Bits")
		ax.set_ylabel("Nr. hits")
		ax.legend()
		ax.set_xlim(xlo, xhi)
		ax.set_ylim(ylo, yhi)
		ax.set_title(title)

ndoms = len(doms)
fig, axs = plt.subplots(nrows=2, ncols=ndoms, figsize=(plotw, ploth))

for domidx in range(ndoms):
	dom = doms[domidx]
	ax = axs[0][domidx]

	xlo = -10
	xhi = 25
	ylo = 0
	yhi = 900
	plot(ax, dom, bin_mids, dom2scop_class2counts[dom], xlo, xhi, ylo, yhi, titles[domidx])

	ax = axs[1][domidx]

	xlo = 0
	xhi = 20
	ylo = 0
	yhi = 100
	plot(ax, dom, bin_mids, dom2scop_class2counts[dom], xlo, xhi, ylo, yhi, None)

plt.tight_layout()
print(plotfn)
plt.savefig(plotfn)
