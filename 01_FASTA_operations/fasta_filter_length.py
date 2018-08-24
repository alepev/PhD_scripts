#!/usr/bin/python

import sys

if len(sys.argv)!=3:
	print
	print 'usage: FASTA FILTER <fasta file> <min LEN>'
	print 'removes fasta sequences with length < min LEN'
	print
else:
	fasta= sys.argv[1]
	L=int(sys.argv[2])
	
	outfile= fasta.split('/')[-1][:-3]+'.len'+str(L)+'.fa'

	PARSE=False
	COUNT=False
	NAME=''
	lc=0
	f=open(outfile,'w')
	for line in open(fasta,'rU').readlines():
			if line[0]=='>':
				if NAME!='' and lc>=LEN:
					f.write(NAME)
					f.write(c)
				NAME=''
				c=''
				lc=0
				PARSE=True
				COUNT=False
				x=line.split('LENGTH=')
				if len(x)>1:
					l= int(x[1].split()[0])
					if l<L:
						PARSE=False
					else:
						f.write(line)
				else:
					NAME=line
					COUNT=True
			elif PARSE:
					if not COUNT:
						f.write(line)
					else:
						c+=line
						lc+=len(line.split()[0])
	if NAME!='' and lc>=LEN:
			f.write(NAME)
			f.write(c)			
	f.close()
	