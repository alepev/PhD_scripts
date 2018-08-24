#!/usr/bin/python

f= open('CLEAN.fa','w')
for line in open('LAST_SELECTION.fa').readlines():
	x=line
	if len(x.split('*'))>1:
		x= line.split('*')[0]+'\n'
	f.write(x)
f.close()
