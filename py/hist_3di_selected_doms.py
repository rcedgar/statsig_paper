import sys

scop_class = None
if len(sys.argv) > 1:
    scop_class = sys.argv[1]
    assert scop_class in [ "a", "b", "c", "d" ]

hitsfn = "../big_hits/3di.scop40"
outprefix = "../hist/3di.scop40."
if not scop_class is None:
    outprefix += scop_class + "."

doms = [ "d1itha_/a.1.1.2", "d1p5dx4/d.129.2.1", "d2g5ra_/b.1.1.1", "d1iyjb4/b.40.4.3" ]
domset = set(doms)
dom2counts = {}

'''
1.47e+03        d2hf2a_/c.108.1.10      d1wl8a1/c.23.16.1
1.38e+03        d2hf2a_/c.108.1.10      d1seza1/c.3.1.2
1.34e+03        d1ng6a_/a.182.1.1       d1ng2a2/b.34.2.1
'''

maxscore = 1000

for dom in doms:
    dom2counts[dom] = [0]*(maxscore+1)

def get_fold(scopid):
    flds = scopid.split('.')
    return flds[0] + '.' + flds[1]

def is_tp(scopid_q, scopid_t):
    fold_q = get_fold(scopid_q)
    fold_t = get_fold(scopid_t)
    return fold_q == fold_t

N = 0
n = 0
for line in open(hitsfn):
    N += 1
    if N%10000 == 0:
        if scop_class is None:
            sys.stdout.write("%d %d\r" % (N, n))
        else:
            sys.stdout.write("%c %d %d\r" % (scop_class, N, n))
    flds = line[:-1].split('\t')
    dom = flds[1]
    if not dom in domset:
        continue
    t = flds[2]
    scopid_q = flds[1].split('/')[1]
    scopid_t = flds[2].split('/')[1]
    if is_tp(scopid_q, scopid_t):
        continue
    scop_class_t = scopid_t[0]
    if not scop_class is None:
        if scop_class_t != scop_class:
            continue
    score = int(float(flds[0]))
    if score > maxscore:
        continue
    n += 1
    dom2counts[dom][score] += 1

for dom in doms:
    counts = dom2counts[dom]
    outfn = outprefix + dom.split('/')[0]
    print("%s %d" % (outfn, sum(counts)))
    f = open(outfn, "w")
    for score in range(maxscore+1):
        f.write("%d\t%d\n" % (score, counts[score]))
    f.close()
