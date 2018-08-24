#!/usr/bin/python

import sys

if len(sys.argv)!=3:
	print 'usage: <fasta file> <string>'
	print 'extracts fasta sequence in which <string> appears in >NAME'
else:
	fasta= sys.argv[1]
	S=sys.argv[2]
	
	out= '.'.join(fasta.split('/')[-1].split('.')[:-1])+'.no_'+S+'.fa'
	
	#parse fasta seqs

	f=open(out,'w')
	APPEND=True
	for line in open(fasta).readlines():
		if line[0]=='>':
			if len(line.split(S))==1:
				f.write(line)
				APPEND=True
			else:
				APPEND=False				
		elif APPEND:
			f.write(line)
	f.close()

	
