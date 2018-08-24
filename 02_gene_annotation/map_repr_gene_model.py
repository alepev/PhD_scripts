#!/usr/bin/python

import sys

models= 'TAIR9_representative_gene_models.txt'

if len(sys.argv)<2:
	print 'usage: <target list (AGI loci)>'
else:
	infile= sys.argv[1]
	outfile=infile[:-4]+'.TAIR9mod.txt'
	
#================================================

#load target list (if any provided)

	list={}
	for line in open(infile,'rU').readlines():
		x= line.split()
		if len(x)!=0 and x[0][0]!='#':
			locus= x[0]
			list[locus]=1
			
	Nlist=len(list.keys())
	print
	print '#gene list (AGI loci) loaded (',Nlist,' genes )'

#================================================

#load representative gene models

	MOD={}
	for line in open(models,'rU').readlines()[3:]:
		x=line.split()
		if len(x)>0:
			locus= x[0].split('.')
			if list.get(locus[0],0)!=0:
				MOD[locus[0]]=locus[1] #locus AGI code + representative model
	
	mods=MOD.keys()
	mods.sort()
	Nloci= len(mods) 
	print		
	print '#representative models loaded (',Nloci,' models )'
	if list!={} and Nloci<Nlist:
		print '***WARNING:',Nlist-Nloci,'provided loci with MISSING MODEL***'
		for locus in list.keys():
			if MOD.get(locus,0)==0:
				print locus
	print
	
#================================================

#output file
	
	f=open(outfile,'w')
	for m in mods:
		f.write(m+'.'+MOD[m]+'\n')
	f.close()







