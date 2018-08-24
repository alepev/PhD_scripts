#! /usr/bin/python

import sys

NTorig=sys.argv[1]
AApred=sys.argv[2]
PRED={}
for line in open(AApred):
	if line[0]=='>':
		p=':'.join(line[1:].split(':')[:-1])
		PRED[p]=1
f=open('.'.join(NTorig.split('.')[:-1])+'.missing.fa','w')
READ=False
for line in open(NTorig):
	if line[0]=='>':
		if PRED.get(line[1:].split()[0],0)==1:
			f.write(line)
			READ=True
		else:
			READ=False
	elif READ:
		f.write(line)
f.close()		
		
	