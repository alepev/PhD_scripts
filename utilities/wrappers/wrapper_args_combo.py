#!/usr/bin/python

import sys,copy
import subprocess as sub

if len(sys.argv)!=3:
	print 'usage: runs script with all combinations of provided arguments'
	print 'wrapper <script> <args> \n(format <args>: "arg1a,arg1b,..|arg2a,arg2b,..|arg3a,arg3b,..")'
	
else:
	script= sys.argv[1]
	args= sys.argv[2].split('|')
	ARGS=[]
	for i in args:
		set=[]
		x=i.split(',')
		for j in x:
			set.append(j)
		ARGS.append(set)

	TEST=ARGS[0]
	l=1
	while l<len(ARGS):
		T2=[]
		for t in TEST:
			for arg in ARGS[l]:
				T2.append(t+' '+arg)
		TEST=copy.copy(T2)
		l+=1
	for args in TEST:
		#print script, args
		print
		sub.call('echo "'+script+' '+args+'";'+script+' '+args,shell=True)
		print
