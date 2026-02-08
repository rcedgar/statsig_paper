import fasta

letters_3di = [ 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' ]
letters_ss = [ 'h', 's', 't', '~' ]

fn_ss = "../data/scop40.ss"
fn_3di = "../data/scop40.3di.fa"

counts = {}
for letter_3di in letters_3di:
	counts[letter_3di] = {}
	for letter_ss in letters_ss:
		counts[letter_3di][letter_ss] = 0

dict_ss = fasta.ReadSeqsDict(fn_ss)
dict_3di = fasta.ReadSeqsDict(fn_3di)

labels = list(dict_ss.keys())
for label in labels:
	seq_ss = dict_ss.get(label)
	seq_3di = dict_3di.get(label)
	if seq_ss is None or seq_3di is None:
		continue
	L = len(seq_ss)
	if len(seq_3di) != L:
		continue
	for i in range(L):
		letter_3di = seq_3di[i].upper()
		letter_ss = seq_ss[i]
		counts[letter_3di][letter_ss] += 1

for letter_3di in letters_3di:
	s = letter_3di
	nh = counts[letter_3di]['h']
	ns = counts[letter_3di]['s']
	nt = counts[letter_3di]['t']
	nl = counts[letter_3di]['~']
	n = nh + ns + nt + nl
	s += "\t%.3g" % (nh/n)
	s += "\t%.3g" % (ns/n)
	s += "\t%.3g" % (nt/n)
	s += "\t%.3g" % (nl/n)
	print(s)

for letter_3di in letters_3di:
	s = "add('" + letter_3di
	nh = counts[letter_3di]['h']
	ns = counts[letter_3di]['s']
	nt = counts[letter_3di]['t']
	nl = counts[letter_3di]['~']
	n = nh + ns + nt + nl
	s += "', fh=%.4f" % (nh/n)
	s += ", fs=%.4f" % (ns/n)
	s += ", ft=%.4f" % (nt/n)
	s += ", fl=%.4f" % (nl/n)
	s += ")"
	print(s)
