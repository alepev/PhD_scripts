#!/usr/bin/python

import sys,copy
import subprocess as sub

if len(sys.argv)!=2:
	print 'usage: runs script with arguments from list in <list_file>'
	print 'wrapper <list_file> \n(list file: one set of arguments per line, args separated by space or tab - use # to exclude lines)'
	
else:
	list_file= sys.argv[1]
	TEST=[]
	for line in open(list_file).readlines():
		x=line.split()
		if len(x)>0 and x[0][0]!='#':
			TEST.append(' '.join(x))
	sub.call('echo "#======================================================================"',shell=True)
	for args in TEST:
		sub.call('echo '+args+'',shell=True)
		sub.call(args,shell=True)
		sub.call('echo "#======================================================================"',shell=True)
