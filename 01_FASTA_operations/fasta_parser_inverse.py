#!/usr/bin/python

import sys

if len(sys.argv)!=3:
	print
	print 'usage: FASTA PARSER <fasta file> <loci list>'
	print 'extracts fasta sequence for loci NOT in list (excludes all models!)'
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
				MOD[x[0]]=1
	
	#load loci list

	LOCI={}
	NOTFOUND=[]
	for line in open(listfile,'rU').readlines():
		x=line.split()[0]
		if len(x)>0 and x[0]!='#':
				if x[0]=='>':
					locus=x.split('>')[1]
				else:
					locus=x
				if locus[:2]=='AT':
					locus=locus.split('|')[0].split('/')[0] #remove descriptors, add here other separators if needed
					LOCI[locus.split('.')[0]]=1
					
	#extract fasta sequence 

	PARSED={}
	READ=False
	for line in open(fasta,'rU').readlines():
			x=line.split()[0]
			if x[0]=='>':
				READ= False
				locus=x[1:]
				if LOCI.get(locus.split('.')[0],0)==0 and MOD.get(locus,0)!=0:
					READ=True
					print line,
			elif READ:
					print line,