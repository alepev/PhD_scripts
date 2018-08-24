#!/usr/bin/python

import sys

if len(sys.argv)!=3:
	print 'usage: <fasta file> <map file>'
	print 'in map file, column1 == ID, column2 == description to add (all other columns ignored)'
	print 'NOTE: no spaces in description! use underscore instead'
else:
	infile=	sys.argv[1]
	mapfile=sys.argv[2]

	M={}
	for line in open(mapfile).readlines():
		if line[0]!='#' and len(line)>2:
			x=line.split()
			if len(x)>=2:
				M[x[0]]=x[1]

	for line in open(infile).readlines():
		if line[0]=='>':
			locus= line[1:].split()[0].split('|')[0]
			if M.get(locus,0)!=0:
				print '>'+locus+'|'+M[locus]
			else:
				print line,
		else:
			print line,

	
