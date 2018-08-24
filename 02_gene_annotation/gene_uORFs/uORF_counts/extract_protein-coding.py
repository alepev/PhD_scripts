#!/usr/bin/python

#WHAT IS THIS PROGRAM FOR:
#extracts protein-coding genes from list

#================================================

#CONTROL PANEL

modelfile=	'TAIR10_representative_gene_models.txt'
gff3=			'TAIR10_GFF3_genes.gff'#'TAIR10_GFF3_genes_transposons.gff'

#================================================

#FUNCTIONS DEFINITION

import sys

#================================================

#START

if len(sys.argv)!=2:
	print
	print 'usage: <target list>'
	print
else:
	infile= sys.argv[1]

#================================================

#load target list + gene models
	
	MOD={}
	for line in open(modelfile,'rU').readlines():
			x=line.split()
			if len(x)!=0 and x[0][:2]=='AT':
				locus=x[0].split('.')
				MOD[locus[0]]=locus[1]
	
	LOCI={} 		
	for line in open(infile,'rU').readlines():
		x=line.split()[0]
		if len(x)>0 and x[0]!='#':
			if x[0]=='>':
				x=x[1:]
			if x[:2]=='AT':
				locus=x.split('|')[0].split('/')[0] #remove descriptors, add here other separators if needed
				if not len(locus.split('.'))>1: #gene model provided
					locus=locus+'.'+MOD[locus]
				LOCI[locus]=1
				
#================================================

#count protein-coding

	PC={}

	for line in open(gff3,'rU').readlines()[1:]:
		x=line.split()
		feat=x[2]
		if feat not in ['five_prime_UTR','three_prime_UTR']:
			locus= x[8].split('=')[1].split(';')[0]
			if LOCI.get(locus,0)!=0:
				if feat=='mRNA':
					PC[locus]=1
	
#================================================
	
#RESULTS

	L= PC.keys()
	L.sort()
	for l in L:
		print l
		 