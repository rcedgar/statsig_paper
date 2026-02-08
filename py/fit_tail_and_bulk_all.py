import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats

SCOP40_DBSIZE = 11211
SCOP40C_DBSIZE = 8340

matplotlib.use("Agg")

plotfn = "../figures/fit_tail_and_bulk_all.svg"
pyfn = "../py/fitted_p_value_params.py"

fpy = open(pyfn, "w")
fpy.write('''import math

name_ref2x1 = {}
name_ref2m0 = {}
name_ref2c0 = {}
name_ref2m = {}
name_ref2c = {}

def set_pvalue_params(name, ref, x1, m0, c0, m, c):
	name_ref2x1[(name, ref)] = x1
	name_ref2m0[(name, ref)] = m0
	name_ref2c0[(name, ref)] = c0
	name_ref2m[(name, ref)] = m
	name_ref2c[(name, ref)] = c

def get_pvalue(name, ref, score):
	x1 = name_ref2x1[(name, ref)]
	if score < x1:
		m = name_ref2m0[(name, ref)]
		c = name_ref2c0[(name, ref)]
	else:
		m = name_ref2m[(name, ref)]
		c = name_ref2c[(name, ref)]
	p = math.exp(m*score + c)
	assert p >= 0
	if p > 1:
		p = 1
	return p

def get_evalue(name, score, dbsize):
	return get_pvalue(name, score)*dbsize

''')

N = SCOP40_DBSIZE
NC = SCOP40C_DBSIZE

NFP_SCOP40 = 125220544
NFP_SCOP40C = 69248284

plotw = 14
ploth = 16

color_counts_fill = "aliceblue"
color_counts_line = "royalblue"
color_fit0 = "darkorchid"
color_fit = "darkorchid"

ls_observed = "dashdot"
ls_fit0 = ":"
ls_fit = "-"

lw_fit0 = 3
lw_fit = 2

lw_counts = 1

titles = [ "FP diff. fold", "FP diff. SF", "FP diff. SF(c)" ]
refs = [ "fold", "sf", "sfc" ]
plotfn = "../figures/fit_tail_and_bulk_all.svg"

hist_dir = "../hist/"
hist_fns = [ "tm.scop40", "foldseek_bits.scop40", "reseek_ts.scop40", "dali.scop40" ]
short_names = [ "tm", "bits", "ts", "dali" ]
axis_labels = [ "TM", "Bits", "ts", "Z" ]
long_names = [ "TM-score", "Foldseek bits", "Reseek ts", "DALI Z" ]

#  0          1  2     3    4        5         6        7
# tm	bin_end	sf	fold	c	cum_sf	cum_fold	cum_c
tsvfldidxs = [ 6, 5, 7 ]

def fit_and_plot(row):
	ax_vec = axs[row]
	assert len(ax_vec) == 3
	short_name = short_names[row]
	long_name = long_names[row]
	axis_label = axis_labels[row]
	hist_fn = hist_dir + hist_fns[row]

	if short_name == "tm":
		#                Fold       SF      SF(c)
		x0s = np.array([ 0.3,		0.3,	0.3], dtype=np.float32)
		y0s = np.array([ 3e-1,		3e-1,	3e-1], dtype=np.float32)

		x1s = np.array([ 0.55,		0.6,	0.6], dtype=np.float32)
		y1s = np.array([ 3e-4,		4e-4,	3e-4], dtype=np.float32)

		x2s = np.array([ 0.78,		0.8,	0.79], dtype=np.float32)
		y2s = np.array([ 1e-7,		5e-7,	1e-7], dtype=np.float32)

		xlos = [ 0.3, 0.3, 0.3 ]
		xhis = [ 0.8, 0.8, 0.8 ]

		ylos = [ 1e-6, 1e-6, 1e-6]
		yhis = [ 1, 1, 1]

	elif short_name == "bits":
		#                Fold       SF      SF(c)
		x0s = np.array([ 20,		20,		20 ], dtype=np.float32)
		y0s = np.array([ 1e-2,		1e-2,	1e-2], dtype=np.float32)

		x1s = np.array([ 80,		100,	100 ], dtype=np.float32)
		y1s = np.array([ 3e-5,		8e-5,	3e-5], dtype=np.float32)

		x2s = np.array([ 350,		350,	170 ], dtype=np.float32)
		y2s = np.array([ 1e-7,		3e-7,	1e-6 ], dtype=np.float32)

		xlos = [20, 20, 20 ]
		xhis = [300, 300, 220 ]

		ylos = [ 1e-6, 1e-6, 1e-6]
		yhis = [ 1e-2, 1e-2, 1e-2]

	elif short_name == "ts":
		#                Fold       SF      SF(c)
		x0s = np.array([ 0.07,		0.07,	0.05], dtype=np.float32)
		y0s = np.array([ 1e-2,		1e-2,	1e-2], dtype=np.float32)

		x1s = np.array([ 0.15,		0.15,	0.11], dtype=np.float32)
		y1s = np.array([ 6e-6,		2e-5,	8e-5], dtype=np.float32)

		x2s = np.array([ 0.35,		0.4,	0.275], dtype=np.float32)
		y2s = np.array([ 1e-8,		1e-8,	1.5e-8], dtype=np.float32)

		xlos = [ 0.025, 0.025, 0.025 ]
		xhis = [ 0.4, 0.4, 0.3 ]

		ylos = [ 1e-8, 1e-8, 1e-8]
		yhis = [ 1e-1, 1e-1, 1e-1]

	elif short_name == "dali":
		x0s = np.array([ 2,		2,		2 ], dtype=np.float32)
		y0s = np.array([ 0.03,	0.03,	0.03], dtype=np.float32)

		x1s = np.array([ 10,		10,		10], dtype=np.float32)
		y1s = np.array([ 3e-5,		3e-4,	2e-4 ], dtype=np.float32)

		x2s = np.array([ 25,		25,		20], dtype=np.float32)
		y2s = np.array([ 1e-7,		1e-6,	5e-7 ], dtype=np.float32)

		xlos = [2,	2, 2 ]
		xhis = [30, 30, 30 ]

		ylos = [ 1e-7, 1e-7, 1e-7]
		yhis = [ 0.1, 0.1, 0.1 ]
	else:
		assert False, f"{short_name=}"

	f = open(hist_fn)
	cmt = f.readline()
	assert cmt.startswith('#')
	hdr = f.readline()
	assert hdr.startswith(short_name)
	cumcounts_vec = []
	scores = []
	cumcounts_vec = []
	for idx in range(3):
		cumcounts_vec.append([])

	for line in f:
		flds = line[:-1].split('\t')
		scores.append(float(flds[0]))
		for idx in range(3):
			cumcounts_vec[idx].append(int(flds[tsvfldidxs[idx]]))

	for idx in range(3):
		ax = ax_vec[idx]
		cumcounts = np.array(cumcounts_vec[idx], dtype=np.int32)
		if idx == 2:
			pvalues = cumcounts/NFP_SCOP40C
		else:
			pvalues = cumcounts/NFP_SCOP40

		x0 = x0s[idx]
		y0 = y0s[idx]

		x1 = x1s[idx]
		y1 = y1s[idx]

		x2 = x2s[idx]
		y2 = y2s[idx]

		log_y0, log_y1 = np.log(y0), np.log(y1)
		m0 = (log_y0 - log_y1) / (x0 - x1)
		c0 = log_y0 - m0 * x0

		log_y1, log_y2 = np.log(y1), np.log(y2)
		m = (log_y2 - log_y1) / (x2 - x1)
		c = log_y1 - m * x1

		x0_vals = np.linspace(x0, x1, 100)
		y0_vals = np.exp(m0 * x0_vals + c0)

		x_vals = np.linspace(x1, x2, 100)
		y_vals = np.exp(m * x_vals + c)

		xlo = xlos[idx]
		xhi = xhis[idx]
		width = (xhi - xlo)/200
	#	ax.bar(scores, pvalues, width=0.005, color=color_counts_fill)
		ax.plot(scores, pvalues, label="Observed", lw=lw_counts, ls=ls_observed, color=color_counts_line)
		ax.fill_between(scores, pvalues, color=color_counts_fill)

		ax.plot(x0_vals, y0_vals, color=color_fit0, ls=ls_fit0, lw=lw_fit0,
			 label=f"Fit ${short_name}<{x1:.2g}$ $p=exp({m0:.2g}*{short_name} {c0:+.2g})$")
		ax.plot(x_vals, y_vals, color=color_fit, ls=ls_fit, lw=lw_fit,
			 label=f"Fit ${short_name}\\geq{x1:.2g}$ $p=exp({m:.2g}*{short_name} {c:+.2g})$")

		fpy.write(f"set_pvalue_params(name='{short_name}', ref='{refs[idx]}', x1={float(x1):.2g}, {m0=:.2g}, {c0=:.2g}, {m=:.2g}, {c=:.2g})\n")

		ax.set_ylabel("p-value")
		ax.set_yscale('log')
		ax.set_xlabel(axis_label)
		ax.set_xlim(xlo, xhi)
		ax.set_ylim(ylos[idx], yhis[idx])
		ax.set_title(long_name + " " + titles[idx])
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3))

fig, axs = plt.subplots(nrows=4, ncols=3)
for row in range(4):
	fit_and_plot(row)
fpy.close()

fig.set_size_inches(w=plotw, h=ploth)
fig.tight_layout()

print(plotfn)
fig.savefig(plotfn)
