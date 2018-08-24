#! /usr/bin/python

mygenes='allgenesArray.txt'
GOann= 'ATH_GO_GOSLIM.txt'

OK={}
for line in open(GOann,'rU').readlines():
	locus= lines.split()[0]
	OK[locus]= OK.get(locus,[]) + [line.split('\t')[5]]

for line in open(mygenes,'rU').readlines()[1:]:
	locus= line.split()[0]
	if OK.get(locus,0)!=0:
		print locus,
		for GO in OK[locus]:
			print '\t'+GO,
		print
