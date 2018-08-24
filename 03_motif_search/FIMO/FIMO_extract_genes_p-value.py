#! /usr/bin/python

import sys

CUT=1000000000

if len(sys.argv)<4 or len(sys.argv)>5: 
	print 'usage: <FIMO output (.TXT)> <test genes> <p-value cutoff> <MAX_DIST (optional)>'
	print
	print 'test presence of motif in provided list and returns genes with p-value for the hit <= given cutoff' 
	#note: selection can be done with q-value, see commented variables in the code
	print 'additionally, a motif localization threshold (i.e. max dist. of motif end from sequence start) can be specified'
	print '(representative gene models are not considered here, since FIMO analysis should be performed already on the representative models)'
	print
else:
	fimofile= sys.argv[1]
	listfile= sys.argv[2]
	P= float(sys.argv[3])
	if len(sys.argv)==5:
		CUT=int(sys.argv[4])
	
	#load target loci list
	
	LOCI={}
	for line in open(listfile,'rU').readlines():
		if len(line)>0 and line[0]!='0':
			locus=line.split('.')[0]
			if line[0]=='>':
				LOCI[locus[1:]]=1
			elif len(locus)==9:
				LOCI[locus]=1
			
	Nloci=len(LOCI.keys())
	
	#count FIMO hits below p-value threshold (and within length limit, if provided)
	
	FOUND={}
	for line in open(fimofile,'rU').readlines()[1:]:	
		x= line.split()
		locus= x[1].split('.')[0]
		start=int(x[2])
		stop=int(x[3])
		pval=float(x[6])
		#qval=float(x[7])
		if stop > CUT or pval > P:
			continue
		elif LOCI.get(locus,0)!=0:
			FOUND[locus]=1
			
	Nfound=len(FOUND.keys())
	Fsort=FOUND.keys()
	Fsort.sort()

	#fold enrichment and Fisher p-value
	
	print
	print '#FIMO file:\t',fimofile
	print '#target loci:\t',listfile
	print '#P-VALUE CUTOFF:\t',P
	print '#DIST.CUTOFF:\t',CUT
	print
	print '#tot loci containing motif:\t',Nfound,'/',Nloci
	print
	for locus in Fsort:
		print locus