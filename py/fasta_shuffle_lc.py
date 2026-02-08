import sys
import fasta
import random

fn = sys.argv[1]
seqs = fasta.ReadSeqsDict(fn)

labels = list(seqs.keys())
nrseqs = len(labels)
seqs = [ seqs[label] for label in labels ]

allseq = "".join(seqs)
allseqv = []

lc = ""
for c in allseq:
    if c.islower():
        lc += c
    else:
        if not lc == "":
            allseqv.append(lc)
        allseqv.append(c)
        lc = ""
if not lc == "":
    v.append(lc)
random.shuffle(allseqv)
alls = "".join(allseqv)

offset = 0
for i in range(nrseqs):
    n = len(seqs[i])
    shuffled_seq = "".join(alls[offset:offset+n])
    assert len(shuffled_seq) == n
    fasta.WriteSeq(sys.stdout, shuffled_seq, labels[i])
    offset += n