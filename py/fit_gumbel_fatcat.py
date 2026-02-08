import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
from lookup import read_lookup

matplotlib.use('Agg')

parser = argparse.ArgumentParser(description="Fit Gumbel distribution with dual-scale visualization.")
parser.add_argument("--input", default="../big_hits/fatcat_len160_q_t_normscore_pvalue")
parser.add_argument("--step", type=int, default=20)
parser.add_argument("--logmaxy", type=float, help="Truncate log y-axis at logmaxy")
parser.add_argument("--logminy", type=float, help="Truncate log y-axis at logminy")
parser.add_argument("--discard_zeros", default=False, action="store_true")
parser.add_argument("--plot", default="../figures/fatcat_gumbel.svg")
parser.add_argument("--hist")
parser.add_argument("--head", type=int)

args = parser.parse_args()

dom2sf, dom2fold = read_lookup("../data/scop40.lookup")

adom = "d1b8da_"

adom_scores_all = []
adom_scores_a = []
adom_scores_b = []

plotw = 10
ploth = 3

scores = []
for line in open(args.input):
	flds = line[:-1].split('\t')
	score = int(float(flds[2]))
	q = flds[0]
	t = flds[1]
	qfold = dom2fold.get(q)
	tfold = dom2fold.get(t)
	if qfold is None or tfold is None:
		continue
	if qfold == tfold:
		continue
	scores.append(score)
	if q == adom:
		adom_scores_all.append(score)
		if tfold.startswith("a."):
			adom_scores_a.append(score)
		elif tfold.startswith("b."):
			adom_scores_b.append(score)

print("%d scores" % len(scores))
print("min %d, max %d" % (np.min(scores), np.max(scores)))

scores = np.array(scores, dtype=np.int32)
mu, beta = stats.gumbel_r.fit(scores)

## params_gevd = stats.genextreme.fit(scores)

mu_fat = 168
beta_fat = 76.1

xlo = 20
xhi = 2000

step = args.step
x_plot_bins = np.arange(xlo, xhi, step)
x_smooth = np.arange(xlo, xhi, step)
y_smooth = stats.gumbel_r.pdf(x_smooth - step/2, loc=mu, scale=beta)
y_smooth_fat = stats.gumbel_r.pdf(x_smooth - step/2, loc=mu_fat, scale=beta_fat)

#  mu_g, beta_g, c_g = params_gevd
mu_g = mu
beta_g = beta
c_g = -0.06
## y_smooth = stats.genextreme.pdf(x_smooth - step/2, loc=mu, scale=beta, c=c_g)

fig, axs = plt.subplots(1, 3, figsize=(plotw, ploth))

for idx in [ 0, 1 ]:
	ax = axs[idx]
	ax.hist(scores, bins=x_plot_bins, width=step*3/4, density=True, color='gray', alpha=0.4, label='Measured')
#	ax.plot(x_smooth, y_smooth, color="black", lw=1, label='EVD fit')
	ax.plot(x_smooth, y_smooth_fat, color="black", lw=1, ls="--", label='FATCAT EVD')
	ax.set_ylabel("PDF")
	ax.set_xlabel("FATCAT score")
	if idx == 0:
		ax.legend()
		ax.set_xlim(0, 800)
		ax.set_title("(a) FPs (fold) $L{\\sim}160$ bulk", fontsize=10)
	if idx == 1:
		ax.set_yscale('log')
		ax.set_xlim(500, 1400)
		ax.set_ylim(1e-9, 1e-3)
		ax.set_title("(b) FPs (fold) $L{\\sim}160$ tail", fontsize=10)

counts_all, _ = np.histogram(adom_scores_all, x_smooth)
counts_a, _ = np.histogram(adom_scores_a, x_smooth)
counts_b, _ = np.histogram(adom_scores_b, x_smooth)

ax = axs[2]
color_all = "gray"
color_a = "green"
color_b = "orange"
alpha = 0.4
xs = x_smooth[:-1]
ax.plot(xs, counts_all, lw=1, color=color_all, alpha=alpha)
ax.fill_between(xs, counts_all, color=color_all, alpha=0.2, label="Total")

ax.plot(xs, counts_a, lw=1, color=color_a, alpha=alpha)
ax.fill_between(xs, counts_a, color=color_a, alpha=0.2, label="All-$\\alpha$")

ax.plot(xs, counts_b, lw=1, color=color_b, alpha=alpha)
ax.fill_between(xs, counts_b, color=color_b, alpha=0.2, label="All-$\\beta$")

ax.set_xlim(50, 800)
ax.set_ylabel("Nr. hits")
ax.set_xlabel("FATCAT score")
ax.set_title("(c) d1b8da_ (all-$\\alpha$)", fontsize=10)
ax.legend()

plt.tight_layout()
plt.savefig(args.plot)
