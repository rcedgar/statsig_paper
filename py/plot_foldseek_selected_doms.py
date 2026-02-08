import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats

matplotlib.use('Agg')

lookupfn = "../data/scop40.lookup"
plotfn = "../figures/foldseek_selected_doms.svg"
plotw, ploth = 9, 4

barw = 3

xlo = -20
xhi = 200
xstep = 5

doms = []
mus = []
lams = []
folds = []

def add_dom(dom, fold, mu, lam):
	doms.append(dom)
	folds.append(fold)
	mus.append(mu)
	lams.append(lam)

add_dom("d1itha_", "a.1", -2.72, 0.156)
add_dom("d1iyjb4", "b.40", -2.37, 0.18)
add_dom("d1p5dx4", "d.129", -2.23, 0.194)

ndoms = len(doms)

dom2sf = {}
for line in open(lookupfn):
	dom, scopid = line[:-1].split('\t')
	flds = scopid.split('.')
	dom2sf[dom] = flds[0] + "." + flds[1] + "." + flds[2]

fig, axs = plt.subplots(nrows=2, ncols=ndoms, figsize=(plotw, ploth))

bin_edges = np.arange(xlo, xhi, xstep)
hist_xs = bin_edges[:-1]

def plot_dom(idx):
	dom = doms[idx]
	mu = mus[idx]
	lam = lams[idx]
	tp_sf = dom2sf[dom]
	fn = "../foldseek_align_selected_doms/" + dom + ".merged"
	f = open(fn)
	line = f.readline()
	assert line.startswith("dom\tdiffscore")
	diffscores = []
	ntp = 0
	for line in f:
		flds = line.split('\t')
		t = flds[0]
		tsf = dom2sf.get(t)
		if tsf is None:
			continue
		if tsf == tp_sf:
			ntp += 1
			continue
		diffscores.append(int(flds[1]))
	print(f"{dom} {tp_sf=} {ntp=} nfp={len(diffscores)}")

	hist, bin_edges2 = np.histogram(diffscores, bin_edges, density=True)
	assert bin_edges2.all() == bin_edges.all()

	beta = 1/lam
	ys = stats.gumbel_r.pdf(hist_xs, loc=mu, scale=beta)

	ax_lin = axs[0][idx]
	ax_log = axs[1][idx]

	for ax in ax_lin, ax_log:
		bar_handle = ax.bar(hist_xs, hist, align="edge", color="darkgray", width=3, label="Observed")
		line_handle = ax.plot(hist_xs, ys, color="black", label="Gumbel", lw=1, ls="dashed")

	ax_lin.set_title(f"{dom} ({folds[idx]})")
	ax_lin.set_ylabel("PDF (linear scale)")
	ax_lin.legend(handles = [bar_handle, line_handle[0]])

	annot_text =  f"$\\mu={mu:.2f}$\n$\\lambda={lam:.2f}$"
	ax_lin.text(0.95, 0.3, annot_text, transform=ax_lin.transAxes, fontsize=10,
		verticalalignment='bottom', horizontalalignment='right')

	ax_log.set_ylabel("PDF (log scale)")
	ax_log.set_yscale('log')
	ax_log.set_ylim(1e-6, 0.1)
	ax_log.set_xlabel('Foldseek "diff." score')

for idx in range(ndoms):
	plot_dom(idx)

plt.tight_layout()
print("Writing " + plotfn)
plt.savefig(plotfn)
