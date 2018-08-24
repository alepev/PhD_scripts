#!/usr/bin/python

import sys

modelfile='TAIR10_representative_gene_models.txt'

if len(sys.argv)!=3:
	print
	print 'usage: <fasta> <CUT>'
	print
else:

	fasta=sys.argv[1]
	CUT=int(sys.argv[2])

	MOD={}
	for line in open(modelfile,'rU').readlines():
			x=line.split()
			if len(x)!=0 and x[0][:2]=='AT':
				MOD[x[0]]=1
	
	NAME=''
	seq=''
	short=1000000000
	long=0
	READ=False
	TOTseqs=0
	BADseqs=0
	for line in open(fasta,'rU').readlines():
			if line[0]=='>':
				if NAME!='' and READ:
						if len(seq)>=CUT:
							if len(seq)>long:
								long=len(seq)
							print '>'+NAME
							print seq[:CUT]
						else:
							BADseqs+=1
							if len(seq)<short:
								short=len(seq)
				NAME=line.split()[0][1:]
				seq=''
				if MOD.get(NAME,0)!=0:
					TOTseqs+=1
					READ=True
				else:
					READ=False
			elif READ:
					seq+=line.split()[0].upper()
	if NAME!='' and READ:
		if len(seq)>=CUT:
			if len(seq)<short:
				short=len(seq)
			print '>'+NAME
			print seq[:CUT]
		else:
			BADseqs+=1
			if len(seq)>long:
				long=len(seq)
				
	#output statistics
	
	print TOTseqs
	print BADseqs
	#print short, long
