#!/usr/bin/python

import sys,copy
import subprocess as sub
sub.PIPE=1

blast= 'blastp'

if len(sys.argv)!=4:
	print
	print 'usage: revblastwrapper.py <e-value> <queries_list> <db (subject)>'
	print '(queries list: one file per line, use # to exclude lines;'
	print 'assumes subject formatted with makeblastdb for -db option)'
	print

else:
	eval= sys.argv[1]
	list_file= sys.argv[2]
	db= sys.argv[3]
	x= sub.Popen('sysctl -n hw.ncpu', shell=True, stdout=1)
	nproc= x.communicate()[0]
	
	FILES=[]
	for line in open(list_file).readlines():
		x=line.split()
		if len(x)>0 and x[0][0]!='#':
			f= x[0]
			FILES.append(f)
			
	for f in FILES:
		out= '.'.join(f.split('.')[:-1])+'.rev'+blast+'.out'
		B= blast+' -db '+db+' -query '+f+' -evalue '+eval+' -out '+out+' -num_threads '+nproc
		sub.call('echo '+B+'',shell=True)
		sub.call(B,shell=True)
		sub.call('echo "#======================================================================"',shell=True)
		