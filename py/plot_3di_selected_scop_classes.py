import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
from viridis_color import viridis_color

## doms = [ "d1itha_/a.1.1.2", "d2g5ra_/b.1.1.1" ]

doms = [ "d1itha_", "d2g5ra_" ]
dom_descs = [ "d1itha_ fold a.1 (all-alpha)", "d2g5ra_ fold b.1 (all-beta)" ]
scop_classes = [ "a", "c", "b" ]
names = [ "all-alpha", "alpha/beta", "all-beta" ]

matplotlib.use('Agg')

plotw = 8
ploth = 5

maxscore = 1000
xlo_lin = 20
xhi_lin = 500
bar_lw = 0
bar_w = 4.5

bin_width = 4
bar_alpha = 0.6

all_color = "gray"

plotfn = "../figures/3di_selected_scop_classes.svg"

scop_class2color = { "a" : viridis_color(0.1), "b" : viridis_color(0.9), "c" : viridis_color(0.5) }

def plot(idx, ax_lin, ax_all):
	dom = doms[idx]
	for classidx, scop_class in enumerate(scop_classes):
		name = names[classidx]
		histfn = "../hist/3di.scop40." + scop_class + "." + dom

		scores_width1 = []
		counts_width1 = []
		with open(histfn) as f:
			for line in f:
				if line.startswith('#'):
					continue
				flds = line[:-1].split('\t')
				scores_width1.append(int(float(flds[0])))
				counts_width1.append(int(float(flds[1])))
		sum_counts = sum(counts_width1)
		scores_width1 = np.array(scores_width1, dtype=np.int32)
		counts_width1 = np.array(counts_width1, dtype=np.int32)
		nbins = len(scores_width1)

		wide_counts = np.array([np.sum(counts_width1[i:i + bin_width]) 
							   for i in range(0, nbins, bin_width)])
		bin_starts = scores_width1[::bin_width]
		bin_mids = bin_starts + bin_width/2

		pdf = wide_counts/(bin_width*sum_counts)

		edgecolor = scop_class2color[scop_class]
		color = scop_class2color[scop_class]
		ax_lin.bar(bin_mids, pdf, linewidth=bar_lw, width=bar_w, alpha=bar_alpha, 
			edgecolor=edgecolor, color=color, label=names[classidx])
		ax_lin.set_xlabel("3Di S-W score")
		ax_lin.set_ylabel("PDF")
		ax_lin.legend()

		ax_lin.set_title(dom_descs[idx])
		ax_lin.set_xlim(xlo_lin, xhi_lin)

	histfn = "../hist/3di.scop40." + dom

	scores_width1 = []
	counts_width1 = []
	with open(histfn) as f:
		for line in f:
			if line.startswith('#'):
				continue
			flds = line[:-1].split('\t')
			scores_width1.append(int(float(flds[0])))
			counts_width1.append(int(float(flds[1])))
	sum_counts = sum(counts_width1)
	scores_width1 = np.array(scores_width1, dtype=np.int32)
	counts_width1 = np.array(counts_width1, dtype=np.int32)
	nbins = len(scores_width1)

	wide_counts = np.array([np.sum(counts_width1[i:i + bin_width]) 
							for i in range(0, nbins, bin_width)])
	bin_starts = scores_width1[::bin_width]
	bin_mids = bin_starts + bin_width/2

	pdf = wide_counts/(bin_width*sum_counts)

	edgecolor = all_color
	color = all_color

	ax_all.bar(bin_mids, pdf, linewidth=bar_lw, width=bar_w, alpha=bar_alpha, 
		edgecolor=edgecolor, color=color, label="all classes")
	ax_all.set_xlabel("3Di S-W score")
	ax_all.set_ylabel("PDF")
	ax_all.legend()

#	ax_all.set_title(dom_descs[idx])
	ax_all.set_xlim(xlo_lin, xhi_lin)

ndoms = len(doms)
fig, axs = plt.subplots(ndoms, 2, figsize=(plotw, ploth))
for domidx in range(ndoms):
	ax_lin = axs[0][domidx]
	ax_all = axs[1][domidx]
	plot(domidx, ax_lin, ax_all)

plt.tight_layout()
print(plotfn)
plt.savefig(plotfn)
