#!/usr/bin/python

import sys,copy
import subprocess as sub
sub.PIPE=1

script= './blast_extract_aligned_region+limits_MERGE_v3_optionalSCORE+FASTA.py'

if len(sys.argv) not in [3,4]:
	print
	print 'usage: wrapper.py <blast_out file list> <extract FASTA? (Y/N)> <SCORE cutoff (optional)>'
	print '(list file: one file per line, use # to exclude lines;'
	print 'assumes subjects formatted with makeblastdb for -db option)'
	print

else:
	list_file= sys.argv[1]
	FASTA= sys.argv[2].upper()
	SCORE= ''
	if len(sys.argv) == 4:
		SCORE= sys.argv[3]
	
	FILES=[]
	for line in open(list_file).readlines():
		x=line.split()
		if len(x)>0 and x[0][0]!='#':
			f= x[0]
			FILES.append(f)
			
	for f in FILES:
		S= script+' '+f+' '+FASTA+' '+SCORE
		print
		sub.call('echo '+S+'',shell=True)
		print
		sub.call(S,shell=True)
		sub.call('echo "#======================================================================"',shell=True)
		