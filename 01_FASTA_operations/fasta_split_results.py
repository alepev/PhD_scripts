#! /usr/bin/python

import sys

FIRST=True
for line in open(sys.argv[1]):
	if line[0]=='>':
		if not FIRST:
			f.close()
		FAM=line[1:4]
		GROUP=line.split('__')[0].split('_')[-1][0]
		f=open(FAM+'_'+GROUP+'.failAug.fa','a')
		f.write(line)
	else:
		f.write(line)
f.close()		
		
	