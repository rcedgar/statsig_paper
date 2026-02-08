import os
import json
import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats

matplotlib.use('Agg')

step = 4
color_gumbel = "black"
color_hist = "gray"
lw_gumbel = 1
plotw = 10
ploth = 5
plotfn_lin = "../figures/3di_lin.svg"
plotfn_log = "../figures/3di_log.svg"
CACHE_DIR = "../tmp_stats_cache"
os.makedirs(CACHE_DIR, exist_ok=True)

def get_cached_stats(fn, scorefld):
	cache_path = os.path.join(CACHE_DIR, os.path.basename(fn) + f".fld{scorefld}.json")
	if os.path.exists(cache_path):
		print(f"Loading cached data for {fn}...")
		with open(cache_path, 'r') as f:
			return json.load(f)
	print(f"Processing raw data for {fn}...")
	df = pd.read_csv(fn, sep='\t', usecols=[scorefld-1], header=None)
	scores = df.iloc[:, 0].values
	mu, beta = stats.gumbel_r.fit(scores)
	counts, bin_edges = np.histogram(scores, bins=np.arange(0, 5000, 4), density=True)
	data = { "mu": float(mu), "beta": float(beta), "counts": counts.tolist(), "bin_edges": bin_edges.tolist() }
	with open(cache_path, 'w') as f:
		json.dump(data, f)
	return data

def do_file(idx, fn, scorefld, log_scale, xlo, xhi, ylo, yhi, title):
	data = get_cached_stats(fn, scorefld)
	mu, beta = data['mu'], data['beta']
	counts = np.array(data['counts'])
	bin_edges = np.array(data['bin_edges'])

	row = idx//3
	col = idx%3
	ax = axs[row][col]

	x_smooth = np.arange(bin_edges[0], xhi, 4)
	y_smooth = stats.gumbel_r.pdf(x_smooth - 2, loc=mu, scale=beta)

	label_gumbel = None
	label_hist = None
	if idx == 0:
		label_gumbel = 'Gumbel fit'
		label_hist = 'Binned scores'

	ax.bar(bin_edges[:-1], counts, width=8, color=color_hist, edgecolor='none',
		alpha=0.4, label=label_hist, align='edge')
	ax.plot(x_smooth, y_smooth, color=color_gumbel, label=label_gumbel, lw=lw_gumbel)
	ax.set_xlabel("3Di S-W score")
	ax.set_ylabel("PDF")

	for patch in ax.patches:
		patch.set_edgecolor('white')
		if log_scale:
			patch.set_linewidth(0.2)
		else:
			patch.set_linewidth(0.3)
	if log_scale:
		annot_text = title
	else:
		annot_text = title + "\n" + r"$\mu$=%.1f $\beta$=%.1f" % (mu, beta)
	ax.text(0.95, 0.95, annot_text, transform=ax.transAxes, fontsize=10,
			verticalalignment='top', horizontalalignment='right')

#	ax.legend(loc='center right', fontsize=10)

	ax.set_xlim(xlo, xhi)
	ax.set_ylim(ylo, yhi)

	print("%d %s mu=%.1f beta=%.1f" % (idx, fn, mu, beta))

	if log_scale:
		ax.set_yscale('log')

plt.rcParams['figure.constrained_layout.use'] = True
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(plotw, ploth))

titles = []
titles.append("Baseline\nAll-vs-all unmodified sequences")
titles.append("FPs (Fold)\nRemoves TP bias\nleaves composition,\nLC and edge bias")
titles.append("Intra-seq. shuffle\nRemoves TP and LC bias\nleaves composition\nand edge bias")
titles.append("Inter-seq. LC-preserve shuffle\nRemoves TP and\ncomposition bias\nleaves LC and\nedge bias")
titles.append("Inter-seq. shuffle\nRemoves TP, LC\nand composition bias\nleaves edge bias")
titles.append("Long inter-seq. shuffle\nRemoves TP, LC,\ncomposition and\nedge bias")

#       idx, fn,                             scorefld, log_scale, xhi, ylo, yhi, title
xlo_lin = 50
xhi_lin = 300
yhi_lin = 0.04
do_file(0, "../big_hits/3di.scop40",				1, False, xlo_lin, xhi_lin, 0, yhi_lin, titles[0])
do_file(1, "../big_hits/3di.scop40.fpsfold",		3, False, xlo_lin, xhi_lin, 0, yhi_lin, titles[1])
do_file(2, "../big_hits/3di.scop40.shuffled",		1, False, xlo_lin, xhi_lin, 0, yhi_lin, titles[2])
do_file(3, "../big_hits/3di.scop40.shuffled_lc",	1, False, xlo_lin, xhi_lin, 0, yhi_lin, titles[3])
do_file(4, "../big_hits/3di.scop40.shuffled_all",	1, False, xlo_lin, xhi_lin, 0, yhi_lin, titles[4])
do_file(5, "../big_hits/3di.scop40.shuffled1k",		1, False, xlo_lin, xhi_lin, 0, yhi_lin, titles[5])

def draw_legend():
	letters=['a', 'b', 'c', 'd', 'e', 'f']
	for idx in range(6):
		row = idx//3
		col = idx%3
		ax = axs[row][col]
		ax.text(0.0, 1.1, letters[idx], transform=ax.transAxes, fontsize=10, fontweight='bold', va='top', ha='left')
	fig.legend(loc='outside lower center', ncols=2)
draw_legend()

#       idx, fn,                             scorefld, log_scale, xhi, ylo, yhi, title

print("Writing " + plotfn_lin)
plt.savefig(plotfn_lin)

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(plotw, ploth))

xlo_log = 0
xhi_log = 1000
ylo_log = 1e-7
yhi_log = 1
do_file(0, "../big_hits/3di.scop40",				1, True, xlo_log, xhi_log, ylo_log, yhi_log, titles[0])
do_file(1, "../big_hits/3di.scop40.fpsfold",		3, True, xlo_log, xhi_log, ylo_log, yhi_log, titles[1])
do_file(2, "../big_hits/3di.scop40.shuffled",		1, True, xlo_log, xhi_log, ylo_log, yhi_log, titles[2])
do_file(3, "../big_hits/3di.scop40.shuffled_lc",	1, True, xlo_log, xhi_log, ylo_log, yhi_log, titles[3])
do_file(4, "../big_hits/3di.scop40.shuffled_all",	1, True, xlo_log, xhi_log, ylo_log, yhi_log, titles[4])
do_file(5, "../big_hits/3di.scop40.shuffled1k",		1, True, xlo_log, xhi_log, ylo_log, yhi_log, titles[5])
draw_legend()

print("Writing " + plotfn_log)
plt.savefig(plotfn_log)
