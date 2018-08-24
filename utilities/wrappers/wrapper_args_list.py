#!/usr/bin/python

import sys,copy
import subprocess as sub

if len(sys.argv)!=3:
	print 'usage: runs script with arguments from list in <list_file>'
	print 'wrapper <list_file> \n(list file: one set of arguments per line, args separated by space or tab - use # to exclude lines)'
	
else:
	script= sys.argv[1]
	list_file= sys.argv[2]
	TEST=[]
	for line in open(list_file).readlines():
		x=line.split()
		set=''
		if len(x)>0 and x[0][0]!='#':
			for arg in x: 
				set+=' '+arg
		if set !='':
			TEST.append(set)
	sub.call('echo "#======================================================================"',shell=True)
	for args in TEST:
		#print script, args
		sub.call('echo '+script+' '+args+'',shell=True)
		sub.call(script+' '+args,shell=True)
		sub.call('echo "#======================================================================"',shell=True)
