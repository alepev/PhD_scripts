#!/usr/bin/python

import sys

for line in open(sys.argv[1],'rU'):
	if line[0]=='>':
		str=line.split()[0]
		L=0
	else:
		if L==0:
			L=len(line)
		else:
			if len(line)!=L:
				print str
