#!/usr/bin/python

import sys

if len(sys.argv)!=3:
	print
	print 'usage: FASTA PARSER <fasta file> <loci list>'
	print 'extracts fasta sequence for loci in list'
	print
else:
	fasta= sys.argv[1]
	listfile=sys.argv[2]
	
	modelfile= 'TAIR10_representative_gene_models.txt'
	outfile= listfile.split('/')[-1].split('.')[0]+'.'+'.'.join(fasta.split('/')[-1].split('.')[:-1])+'.fa'
	
	#load representative gene models
		
	MOD={}
	for line in open(modelfile,'rU').readlines():
			x=line.split()
			if len(x)!=0 and x[0][:2]=='AT':
				locus=x[0].split('.')
				MOD[locus[0]]=locus[1]
	
	#load loci list

	LOCI={}
	NOTFOUND=[]
	for line in open(listfile,'rU').readlines():
		x=line.split()[0]
		if len(x)>0 and x[0]!='#' and x[:2]=='AT':
			locus=x.split('|')[0].split('/')[0] #remove descriptors, add here other separators if needed
			if len(locus.split('.'))>1: #gene model provided
				LOCI[locus]=1
			elif MOD.get(locus,0)!=0: #representative gene model automatically assigned
				LOCI[locus+'.'+MOD[locus]]=1
			else:
				NOTFOUND.append(locus)
	
	if NOTFOUND!=[]:
		print
		print '***WARNING: the following IDs could not be mapped to a gene model (excluded from parsing):'
		for locus in NOTFOUND:
			print locus
		print
					
	#extract fasta sequence 

	PARSED={}
	READ=False
	f=open(outfile,'w')
	for line in open(fasta,'rU').readlines():
			x=line.split()[0]
			if x[0]=='>':
				READ= False
				locus=x[1:]
				if LOCI.get(locus,0)!=0:
					READ=True
					PARSED[locus]=1
					f.write(line)
			elif READ:
					f.write(line)
	f.close()
	
	if len(PARSED)!=len(LOCI):
		print
		print '***WARNING: the following IDs do not have an associated sequence (excluded from output):'
		l=LOCI.keys()
		l.sort()
		for locus in l:
			if PARSED.get(locus,0)==0:
				print locus
		print
	