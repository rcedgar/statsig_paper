import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import gumbel_r
from scipy import stats

# Xu & Zhang parameters from paper
mu_xz = 0.1512
beta_xz = 0.0242

matplotlib.use('Agg')

# Integral of Gumbel distribution from x_start to 1.
def calculate_gumbel_integral(mu, beta, x_start):
	dist = gumbel_r(loc=mu, scale=beta)
	cdf_upper = dist.cdf(1)
	cdf_lower = dist.cdf(x_start)
	return cdf_upper - cdf_lower

def fmt(value):
    s = "{:.2e}".format(value)
    mantissa, exponent = s.split('e')
    exponent = int(exponent) 
    return fr"${mantissa} \times 10^{{{exponent}}}$"

plotfns = [ "../figures/tm_gumbel.png", "../figures/tm_gumbel.svg" ]
plotw = 10
ploth = 3
nrbins = 100
scorefldnr = 3
xlo = 0
xhi = 0.8

ylo_lin = 0
ylo_log = 1e-4

yhi_lin = None
yhi_log = 40

color_same = "lightgray"
color_diff = "darkgray"

fig, axs = plt.subplots(1, 2, figsize=(plotw, ploth))
ax_lin, ax_log = axs

df_all = pd.read_csv("../big_data/tm.scop40", sep='\t', usecols=[scorefldnr-1], header=None)
df_fpfold = pd.read_csv("../big_hits/tm.scop40.fpsfold", sep='\t', usecols=[scorefldnr-1], header=None)

scores_all = df_all.iloc[:, 0].values
scores_fpfold = df_fpfold.iloc[:, 0].values

x_plot = np.linspace(xlo, xhi, 1000)

mu_scipy, beta_scipy = stats.gumbel_r.fit(scores_fpfold)
print(f'{mu_scipy=:.3e}, {beta_scipy=:.3e}')
# mu_scipy=2.525e-01, beta_scipy=4.398e-02

y_xz = stats.gumbel_r.pdf(x_plot, loc=mu_xz, scale=beta_xz)
y_scipy = stats.gumbel_r.pdf(x_plot, loc=mu_scipy, scale=beta_scipy)

ax_lin.plot(x_plot, y_xz, lw=1, color="black", linestyle="solid", label="Xu & Zhang fit")
ax_log.plot(x_plot, y_xz, lw=1, color="black", linestyle="solid", label="Xu & Zhang fit")

ax_lin.plot(x_plot, y_scipy, lw=1, color="black", linestyle="dashed", label="scipy fit")
ax_log.plot(x_plot, y_scipy, lw=1, color="black", linestyle="dashed", label="scipy fit")

ax_lin.hist(scores_all, bins=nrbins, density=True, color=color_same, edgecolor=color_same, label='same fold')
ax_log.hist(scores_all, bins=nrbins, density=True, color=color_same, edgecolor=color_same, label='same fold')

ax_lin.hist(scores_fpfold, bins=nrbins, density=True, color=color_diff, edgecolor=color_diff, label='different fold')
ax_log.hist(scores_fpfold, bins=nrbins, density=True, color=color_diff, edgecolor=color_diff, label='different fold')

ax_lin.set_xlim(xlo, xhi)
ax_log.set_xlim(xlo, xhi)

ax_lin.set_ylim(ylo_lin, yhi_lin)
ax_log.set_ylim(ylo_log, yhi_log)

ax_lin.set_xlabel("TM-score")
ax_log.set_xlabel("TM-score")

ax_lin.set_ylabel("PDF")
ax_log.set_ylabel("PDF")

ax_lin.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax_log.xaxis.set_major_locator(ticker.MultipleLocator(0.1))

pvalue_scipy = calculate_gumbel_integral(mu_scipy, beta_scipy, 0.5)
pvalue_xz = calculate_gumbel_integral(mu_xz, beta_xz, 0.5)
pvalue_measured = np.mean(scores_fpfold >= 0.5)

annot_text_lin = r"Xu & Zhang $\mu$=%.3f $\beta$=%.3f" "\n" r"scipy $\mu$=%.3f $\beta$=%.3f" % \
	(mu_xz, beta_xz, mu_scipy, beta_scipy)

annot_text_log = "$p$-value at TM=0.5\nXu & Zhang fit $p$=%s\nmeasured $p$=%s\nscipy fit $p$=%s" % \
	(fmt(pvalue_xz), fmt(pvalue_measured), fmt(pvalue_scipy))

ax_lin.text(0.95, 0.95, annot_text_lin, transform=ax_lin.transAxes, fontsize=10,
		verticalalignment='top', horizontalalignment='right')

ax_log.text(0.95, 0.95, annot_text_log, transform=ax_log.transAxes, fontsize=10,
		verticalalignment='top', horizontalalignment='right')

ax_lin.legend(loc='lower right')

ax_log.set_yscale('log')

for plotfn in plotfns:
	plt.tight_layout()
	sys.stderr.write(plotfn + '\n')
	plt.savefig(plotfn)
