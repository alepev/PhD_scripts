#!/usr/bin/python

import sys
import subprocess as sub

if len(sys.argv)!=3:
	print 'usage: runs script with arguments from list in <list_file>'
	print 'wrapper <script> <list_file> \n(list file: one set of arguments per line, args separated by space or tab - use # to exclude lines)'
	
else:
	script= sys.argv[1]
	if script[:2]!='./':
		script='./'+script
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
	if TEST!=[]:
		print
		print '# SCRIPT:\t',script
		print
		print '# LIST FILE:\t',list_file
		sub.call('cat '+list_file+'',shell=True)
		print
		print
		for i in range(len(TEST)):
			args=TEST[i]
			print '#======================================================================'
			print 'now running:\t',script,args
			#sub.call('echo '+script+' "'+args+'"',shell=True) #diagnostic
			sub.call(script+' '+args,shell=True)
		print
		print '***DONE***'
	else:
		print '***WARNING: EMPTY LIST PROVIDED!***'
	print