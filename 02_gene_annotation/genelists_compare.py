#!/usr/bin/python

import sys

if len(sys.argv)<4 or len(sys.argv)>5 or sys.argv[3] not in ['U1','U2','OV','ALL']:
	print 'usage: list1, list2, U1/U2/OV/ALL, <out_name> - only if "ALL"'
	print '("U1" or "U2" == unique for list 1 or 2, "OV" == overlap, "ALL" == all of the 3 in separate files)'
else:	
	f1=str(sys.argv[1])
	f2=str(sys.argv[2])
	ov= str(sys.argv[3])
	if ov=='ALL':
		out= str(sys.argv[4])
	
	l1=[]
	for line in open(f1).readlines():
		x=line.split()
		if len(x) > 0:
			l1.append(x[0])
	if ov!= 'ALL':
		if ov=='U1':
			for line in open(f2).readlines():
				x=line.split()
				if len(x) > 0 and x[0] in l1:
					l1.pop(l1.index(x[0]))
		elif ov=='U2':
			l2=[]
			for line in open(f2).readlines():
				x=line.split()
				if len(x) > 0 and x[0] not in l1:
					l2.append(x[0])
			l1=l2
		elif ov=='OV':
			l2=[]
			for line in open(f2).readlines():
				x=line.split()
				if len(x) > 0 and x[0] in l1:
					l2.append(x[0])
			l1=l2
		for x in l1:
			print x
	else:
		OV=[]
		l2=[]
		for line in open(f2).readlines():
			x=line.split()
			if x[0] in l1:
				l1.pop(l1.index(x[0]))
				OV.append(x[0])
			else:
				l2.append(x[0])	

		f=open(out+'_UNIQ1.txt','w')
		for x in l1:
			f.write(x+'\n')
		f.close()	
		f=open(out+'_UNIQ2.txt','w')
		for x in l2:
			f.write(x+'\n')
		f.close()
		f=open(out+'_OVERLAP.txt','w')
		for x in OV:
			f.write(x+'\n')
		f.close()	
