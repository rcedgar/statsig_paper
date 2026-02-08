#!/usr/bin/python3
import sys
import argparse
import numpy as np

# sys.stderr.write(' ' .join(sys.argv) + '\n')
AP = argparse.ArgumentParser()
AP.add_argument("--hits", required=True)
AP.add_argument("--lookup", default="../data/scop40.lookup")
AP.add_argument("--lookupc", default="../data/scop40c.lookup")
AP.add_argument("--fields", type=int, required=True, nargs=3)
AP.add_argument("--minlength", type=int)
AP.add_argument("--maxlength", type=int)
AP.add_argument("--type", type=str, choices=["bits", "tm", "ts", "blastp", "dali"], required=True)
AP.add_argument("--output", type=str, required=True)
args = AP.parse_args()

qfld, tfld, scorefld = args.fields

if args.type == "bits":
	minscore = 20
	maxscore = 5000
elif args.type == "ts":
	minscore = -0.25
	maxscore = 1
elif args.type == "tm":
	minscore = 0.0
	maxscore = 1
elif args.type == "blastp":
	minscore = 40
	maxscore = 4000
elif args.type == "dali":
	minscore = 2
	maxscore = 30

dom2sf = {}
dom2fold = {}
scop40c_doms = set()
for line in open(args.lookupc):
	scop40c_doms.add(line.split('\t')[0].split('/')[0])

for line in open(args.lookup):
	dom, fam = line[:-1].split('\t')
	flds = fam.split('.')
	fold = flds[0] + '.' + flds[1]
	dom2fold[dom] = fold
	dom2sf[dom] = fold + '.' + flds[2]

N = None
Nc = None
if not args.minlength is None:
	N = 0
	Nc = 0
	dom2length = {}
	for line in open("../data/scop40.seqlengths"):
		flds = line[:-1].split('\t')
		assert len(flds) == 2
		dom = flds[0].split('/')[0]
		length = int(flds[1])
		if length >= args.minlength and length <= args.maxlength:
			N += 1
			if dom in scop40c_doms:
				Nc += 1
		dom2length[dom] = length
	sys.stderr.write("N=%d Nc=%d\n" % (N, Nc))

all_scores = []
scores_fpsf = []
scores_fpfold = []
scores_fpc = []
scores_tpc = []
n = 0
for line in open(args.hits):
	flds = line[:-1].split()
	q = flds[qfld-1].split('/')[0]
	t = flds[tfld-1].split('/')[0]
	if q == t:
		continue
	score = float(flds[scorefld-1])
	if score < minscore or score > maxscore:
		continue
	if not args.minlength is None:
		ql = dom2length.get(q)
		tl = dom2length.get(t)
		if ql is None or tl is None:
			continue
		if ql < args.minlength or ql > args.maxlength:
			continue
		if tl < args.minlength or tl > args.maxlength:
			continue
	qsf = dom2sf.get(q)
	tsf = dom2sf.get(t)
	if qsf is None or tsf is None:
		continue
	qfold = dom2fold[q]
	tfold = dom2fold[t]
	in_scop40c = q in scop40c_doms and t in scop40c_doms

	all_scores.append(score)
	if in_scop40c and qsf == tsf:
		scores_tpc.append(score)

	if qfold != tfold:
		scores_fpsf.append(score)
		scores_fpfold.append(score)
		if in_scop40c:
			scores_fpc.append(score)
	elif qsf != tsf:
		scores_fpsf.append(score)
		if in_scop40c:
			scores_fpc.append(score)

scores_fpsf = np.array(scores_fpsf, dtype=np.float32)
scores_fpfold = np.array(scores_fpfold, dtype=np.float32)
scores_fpc = np.array(scores_fpc, dtype=np.float32)

d_min = min(np.min(scores_fpsf), np.min(scores_fpfold), np.min(scores_fpc))
d_max = max(np.max(scores_fpsf), np.max(scores_fpfold), np.max(scores_fpc))

if args.type == "bits" or args.type == "blastp":
	bins = np.arange(d_min, d_max + 2)
elif args.type == "tm" or args.type == "ts" or args.type == "dali":
	bins = np.linspace(minscore, maxscore, 101)
else:
	assert False

all_counts, bin_edges = np.histogram(all_scores, bins=bins)

counts_tpc, _ = np.histogram(scores_tpc, bins=bin_edges)
counts_fpsf, _ = np.histogram(scores_fpsf, bins=bin_edges)
counts_fpfold, _ = np.histogram(scores_fpfold, bins=bin_edges)
counts_fpc, _ = np.histogram(scores_fpc, bins=bin_edges)

n = len(counts_fpsf)
assert len(counts_fpfold) == n
assert len(counts_fpc) == n
for i in range(n):
	if counts_fpfold[i] > counts_fpsf[i]:
		sys.stderr.write("bin %d score %.3g count_fpfold=%d count_fpsf=%d\n" %
				   (i, bin_edges[i], float(bin_edges[i]), counts_fpfold[i], counts_fpsf[i]))

cumsum_fpsf = np.cumsum(counts_fpsf[::-1])[::-1]
cumsum_fpfold = np.cumsum(counts_fpfold[::-1])[::-1]
cumsum_fpc = np.cumsum(counts_fpc[::-1])[::-1]
cumsum_tpc = np.cumsum(counts_tpc[::-1])[::-1]

with open(args.output, "w") as f:
	f.write("# " + ' ' .join(sys.argv))
	if not N is None:
		f.write(" # N=%d" % N)
	if not Nc is None:
		f.write(" # Nc=%d" % Nc)
	f.write("\n")
	f.write("%s\tbin_end\tsf\tfold\tc\tcum_sf\tcum_fold\tcum_c\ttpc\tcum_tpc\n" % args.type)
	for i in range(len(counts_fpsf)):
		if args.type == "bits" or args.type == "blastp":
			s = "%d" % bin_edges[i]
			s += "\t%d" % bin_edges[i+1]
		elif args.type == "tm" or args.type == "ts" or args.type == "dali":
			s = "%.4g" % bin_edges[i]
			s += "\t%.4g" % bin_edges[i+1]
		else:
			assert False
		s += "\t%d" % counts_fpsf[i]
		s += "\t%d" % counts_fpfold[i]
		s += "\t%d" % counts_fpc[i]

		s += "\t%d" % cumsum_fpsf[i]
		s += "\t%d" % cumsum_fpfold[i]
		s += "\t%d" % cumsum_fpc[i]

		s += "\t%d" % counts_tpc[i]
		s += "\t%d" % cumsum_tpc[i]
		f.write(s + '\n')
