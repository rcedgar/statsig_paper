import math
import math
import random

bfn = "../big_hits/blastp.scop40"
rfn = "../big_hits/reseek27_fast.scop40"

sample_size = 1000

'''
# head blastp.scop40
d3nfka_ d3nfka_ 1.14e-65        481
d3nfka_ d1qava_ 2.79e-11        124
d3nfka_ d3i4wa_ 5.60e-07        96
'''
qtset = set()
qt2blastpe = {}
for line in open(bfn):
	flds = line.split()
	assert len(flds) == 4
	q, t = flds[0], flds[1]
	if q == t:
		continue
	e = float(flds[2])
	qtset.add((q,t))
	qt2blastpe[(q,t)] = e

'''
# head reseek27_fast.scop40
d1tdja3/d.58.18.2       d1tdja3/d.58.18.2       2.3e-09
d3nfka_/b.36.1.1        d3nfka_/b.36.1.1        7.3e-11
d1neia_/d.253.1.1       d1neia_/d.253.1.1       5.53e-07
'''
v = []
for line in open(rfn):
	flds = line.split()
	assert len(flds) == 3
	q, t = flds[0], flds[1]
	q = q.split('/')[0]
	t = t.split('/')[0]
	if q == t:
		continue
	if not (q,t) in qtset:
		continue
	be = qt2blastpe[(q,t)]
	re = float(flds[2])
	v.append((q,t,be,re))

random.shuffle(v)
v = v[:sample_size]
for q, t, be, re in v:
	s = q
	s += "\t" +  t
	s += "\t%.3g" % be
	s += "\t%.3g" % re
	print(s)
