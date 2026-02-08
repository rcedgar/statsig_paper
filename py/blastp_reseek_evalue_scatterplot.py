import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from viridis_color import viridis_color

matplotlib.use("Agg")

plotfn = '../figures/blastp_reseek_evalue_scatterplot.svg'
plotw = 4
ploth = 4

bes = []
res = []
for line in open("../results/blastp_reseek_evalues_sample.tsv"):
    flds = line.split()
    assert len(flds) == 4
    be = float(flds[2])
    re = float(flds[3])
    bes.append(be)
    res.append(re)

logbes = np.log10(bes)
logres = np.log10(res)

a, b = np.polyfit(logres, logbes, 1)
xg = np.logspace(np.log10(min(bes)), np.log10(max(bes)), 200)
yg = 10**b * xg**a

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(plotw, ploth))

ax.plot(res, bes, '.', markersize=1, color="gray")
ax.plot(xg, yg, lw=1, color="black", label=f"Log-linear fit")

ax.set_xlabel('Reseek E-value')
ax.set_ylabel('BLASTP E-value')
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend(loc='upper left')
ax.set_xlim(1e-30, 1)
ax.set_ylim(1e-30, 1)
ax.set_aspect('equal', adjustable='box')

plt.tight_layout()
print(plotfn)
plt.savefig(plotfn)
