#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

plotw = 5
ploth = 3
xlo = -300
xhi = 300
step = 4

bin_edges = np.arange(xlo, xhi, step)

hitsfn = "../big_hits/3di.scop40.shuffled_lc_subset"
plotfn = "../figures/fit_logistic.svg"

vals = []
with open(hitsfn, "r", encoding="utf-8", errors="replace") as f:
	for line in f:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		fields = line.split("\t")
		if not fields:
			continue
		try:
			vals.append(float(fields[0]))
		except ValueError:
			continue

scores = np.asarray(vals, dtype=float)
indices_i = np.random.randint(0, len(scores), size=len(scores))
indices_j = np.random.randint(0, len(scores), size=len(scores))
x = scores[indices_i] - scores[indices_j]

loc, scale = stats.logistic.fit(x)  # MLE

xaxis = bin_edges[:-1]
pdf = stats.logistic.pdf(xaxis, loc=loc, scale=scale)

fig = plt.figure(figsize=(plotw, ploth))
ax = plt.gca()

h, _ = np.histogram(x, bins=bin_edges, density=True)
print(f"{len(h)=} {len(xaxis)=} {len(pdf)=} {len(bin_edges)=}")

# ax.hist(x-0.5, bins=bin_edges, density=True, label="Measured", color="lightgray", rwidth=2)
ax.fill_between(xaxis, h, label="Measured", color="lightgray")
ax.plot(xaxis, pdf, linewidth=1, label="Logistic fit", color="black")

ax.set_xlim(-400, 400)
ax.set_xlabel("Random (score$_i$ $-$ score$_j$)")
ax.set_ylabel("PDF")

etext = f"loc={loc:.3g}\nscale={scale:.3g}"
ax.text(150, 0.003, etext, fontsize=10)

#ax.set_title(f"Logistic fit to difference of S-W scores\nloc={loc:.3g}, scale={scale:.3g}")
ax.legend()

fig.tight_layout()
print("Writing " + plotfn)
fig.savefig(plotfn)
plt.close(fig)
