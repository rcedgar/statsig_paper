import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.stats
from scipy import stats

fold_or_sf = "fold"
assert fold_or_sf == "fold" or fold_or_sf == "sf"

matplotlib.use('Agg')

def fmt(value):
    s = "{:.2e}".format(value)
    mantissa, exponent = s.split('e')
    exponent = int(exponent) 
    return fr"${mantissa} \times 10^{{{exponent}}}$"

plotfns = [ "../figures/tm_gevd_figure_bad_scipy_fit.svg" ]

plotw = 10
ploth = 3
nrbins = 100
scorefldnr = 3
xlo = 0
xhi = 0.8

ylo_lin = 0
ylo_log = 1e-2

yhi_lin = None
yhi_log = 40

rng = np.random.default_rng()
np.random.seed(1)

fig, axs = plt.subplots(1, 2, figsize=(plotw, ploth))
ax_lin, ax_log = axs

print("../big_hits/tm.scop40.fps" + fold_or_sf)
df_fpfold = pd.read_csv("../big_hits/tm.scop40.fps" + fold_or_sf, sep='\t', usecols=[scorefldnr-1], header=None)

scores_fpfold = df_fpfold.iloc[:, 0].values

# Fails to converge with full dataset, try subset
scores_fpfold = rng.choice(scores_fpfold, size=10000, replace=False)

x_plot = np.linspace(xlo, xhi, 1000)

mu_manual = 2.525e-01
beta_manual = 4.398e-02
c_manual = 0.05
xi_manual = -c_manual

print("scipy fitting")

params_gevd = stats.genextreme.fit(scores_fpfold)

mu_scipy, beta_scipy , c_scipy = params_gevd
xi_scipy = -c_scipy

print(f'{mu_scipy=:.3g} {beta_scipy=:.3g} {xi_scipy=:.3g}')

y_scipy = stats.genextreme.pdf(x_plot, loc=mu_scipy, scale=beta_scipy, c=-xi_scipy)
y_manual = stats.genextreme.pdf(x_plot, loc=mu_manual, scale=beta_manual, c=-xi_manual)

ax_lin.plot(x_plot, y_scipy, lw=1, color="black", linestyle="dashed", label="Scipy fit")
ax_log.plot(x_plot, y_scipy, lw=1, color="black", linestyle="dashed", label="Scipy fit")

ax_lin.plot(x_plot, y_manual, lw=1, color="black", linestyle="solid", label="Manual fit")
ax_log.plot(x_plot, y_manual, lw=1, color="black", linestyle="solid", label="Manual fit")

color = "lightgray"
ax_lin.hist(scores_fpfold, bins=nrbins, density=True, color=color, edgecolor=color, label='FPs (fold)')
ax_log.hist(scores_fpfold, bins=nrbins, density=True, color=color, edgecolor=color, label='FPs (fold)')

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

ax_log.set_yscale('log')

ax_lin.legend()
ax_log.legend()

for plotfn in plotfns:
	plt.tight_layout()
	sys.stderr.write(plotfn + '\n')
	plt.savefig(plotfn)
