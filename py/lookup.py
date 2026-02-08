def read_lookup(lookupfn):
	dom2sf = {}
	dom2fold = {}
	for line in open(lookupfn):
		dom, scopid = line[:-1].split('\t')
		flds = scopid.split('.')
		sf = flds[0] + "." + flds[1] + "." + flds[2]
		dom2sf[dom] = sf
		fold = flds[0] + "." + flds[1]
		dom2fold[dom] = fold
	return dom2sf, dom2fold
