#! /usr/bin/python

#version 2.0 with % hits over total!

import sys
import os
import copy
import subprocess as sub
sub.PIPE=1

########################################################

def FoldFisher(N,K,n,k):
	#N == tot. genes in background
	#K == tot. genes in background with annotation X
	#n == genes in list
	#k == genes in list with annotation X
	if k==0 or n==0 or K==0 or N==0:
		x,pval= '-','-'
	else:
		FC=(float(k)/float(n))/(float(K)/float(N))
		if FC>=1:
			alt='greater'
		else:
			alt='less'
		y=sub.Popen('Rscript Fisher_1-tailed.R '+str(N)+' '+str(K)+' '+str(n)+' '+str(k)+' '+alt, shell=True, stdout=1)
		x= '%1.2f' % (FC)
		pval='%1.2e' % float(y.communicate()[0])
	return (x,pval)

#ex. of Fisher test for motif m:
#FISHER[m]= FoldFisher(Nbkg,Mbkg[m],Nloci,Mloci[m])

########################################################

if len(sys.argv)<4 or len(sys.argv)>5:
	print 'usage: <FIMO output (.TXT)> <background.FASTA> <target_loci.FASTA> <MAX_DIST (optional)>'
	print
	print 'computes fold enrich. + p-value (1-tailed Fisher exact test) of target list vs background,' 
	print 'using different motif significance thresholds (see code to select threshold);' #you can set thresholds to test modifying "Pthr" below
	print 'additionally, a motif localization requirement (i.e. max distance of motif end from sequence start) can be specified'
	print 
	print 'NOTE: make sure to specify exact background used as FIMO input, and a corresponding target loci file (FASTA format);'
	print 'for instance, if FIMO analysis was done on CDS sequences, provide CDS sequence files for both background and targets.'
	print '(no need to provide representative gene models here, since FIMO analysis should have been performed already for the representative models only)'
	print
else:
	Pthr=[1e-3,1e-4,1e-5,1e-6]
	CUT=1000000000
	fimofile= sys.argv[1]
	bkgfile= sys.argv[2]
	listfile= sys.argv[3]
	if len(sys.argv)==5:
		CUT=int(sys.argv[4])
	
	#load target loci list
	
	LOCI={}
	for line in open(listfile,'rU').readlines():
		if len(line)>0 and line[0]=='>':
			locus= line.split('.')[0][1:]
			LOCI[locus]=1
	
	Nloci=len(LOCI.keys())
	
	#count tot genes in background
	
	x=sub.Popen('cat '+bkgfile+' | grep ">" | wc -l', shell=True, stdout=1)
	Nbkg= int(x.communicate()[0])
	
	#count motif matches
	
	Mbkg,Mloci={},{}
	FOUND={}
	for p in Pthr:
		Mbkg[p],Mloci[p]=0,0
		FOUND[p]={}
	for line in open(fimofile,'rU').readlines()[1:]:	
		x= line.split()
		locus= x[1].split('.')[0]
		start=int(x[2])
		stop=int(x[3])
		pval=float(x[6])
		if stop > CUT:
			continue
		for p in Pthr:
			if pval <= p:
				if FOUND[p].get(locus,0)==0:
					if LOCI.get(locus,0)!=0:
						Mloci[p]=Mloci[p]+1
					else:
						Mbkg[p]=Mbkg[p]+1
					FOUND[p][locus]=1
			else:
				break

	#fold enrichment and Fisher p-value
	
	print
	print '#FIMO file:\t',fimofile
	print '#target loci:\t',listfile
	print '#DIST.CUTOFF:\t',CUT
	for p in Pthr:
		FISHER= FoldFisher(Nbkg,Mbkg[p],Nloci,Mloci[p])
		print
		print '#MOTIF P-VALUE THRESHOLD:\t',p
		print '#TOTloci\tTOTbkg\tNloci\tNbkg\tN%loci\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)'
		print str(Nloci)+'\t'+str(Nbkg)+'\t'+str(Mloci[p])+'\t'+str(Mbkg[p]),
		print '\t%4.2f\t%4.2f' % (100.0*Mloci[p]/Nloci, 100.0*Mbkg[p]/Nbkg),
		print '\t'+FISHER[0]+'\t'+FISHER[1]
	print