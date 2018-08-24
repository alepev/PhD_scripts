#!/usr/bin/python

import sys

for line in open(sys.argv[1]).readlines():
	if line[0]=='>':
		x=line.split()[0]
		if len(x.split('|'))>1:
			x='|'.join(x.split('|')[:-1])
		print x
	else:
		print line,
