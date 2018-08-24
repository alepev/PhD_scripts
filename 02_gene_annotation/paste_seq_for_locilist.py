#!/usr/bin/python

import os,sys

if len(sys.argv)!=4:
	print 'parse sequence for lists in <lists folder> from fasta file <sequence file>, outfile description <suffix>'
	print '(suffix will be added as "_<suffix>.txt" - ex. 500upstream, intronic, etc.'
else:
	targetdir= sys.argv[1]
	seqfile= sys.argv[2]
	if targetdir[-1]!='/':
		targetdir=targetdir+'/'
	suffix=sys.argv[3]
	
#####################################################################
# LOAD upstream sequences
#####################################################################

	seq={}
	for line in open(seqfile).readlines():
		x= line.split()
		if line[0]=='>':
			locus=x[0][1:]
			seq[locus]=''
		else:
			seq[locus]= seq[locus]+str(x[0])

#####################################################################
# LOAD target list and OUTPUT target sequences
#####################################################################

	for t in os.listdir(targetdir):
		outfile= t[:-4]+'_'+suffix+'.txt'
		f=open(outfile,'w')
		c=0
		for line in open(targetdir+t).readlines():
			x= line.split()
			#empty lines ignored + target diff. from over-expr. bZIPs + sequence presence check
			if len(x)>0 and seq.get(x[0],0)!=0: 
				f.write('>'+x[0]+'\n'+seq[x[0]]+'\n')
				c+=1
		f.close()
		if c>0:
			print 'sequence(s) parsed for',t,':\t',c
		else:
			print '***warning: no output sequence parsed for ',t,'!'
			os.remove(outfile)
