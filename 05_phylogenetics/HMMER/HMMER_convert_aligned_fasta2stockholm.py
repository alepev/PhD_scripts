#!/usr/bin/python

import sys

if len(sys.argv)!=2:
	print 'usage: converts <aligned FASTA> in STOCKHOLM format'
else:
	fasta= sys.argv[1]
	stock= '.'.join(fasta.split('.')[:-1])+'.sto'
	
	f= open(stock,'w')
	f.write('# STOCKHOLM 1.0')
	for line in open(fasta).readlines():
		if line.split()!=[]:
			if line[0]=='>':
				f.write('\n'+'_'.join(line[1:-1].split())+' ')
			else:
				f.write(line[:-1].upper())
	f.write('\n//')
	f.close()
