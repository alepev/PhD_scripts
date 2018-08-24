#!/usr/bin/python

import sys
import os
import subprocess as sub
sub.PIPE=1

if len(sys.argv)!=2:
	print
	print 'usage: <file OR folder>'	
	print 'maps PAC_IDs to gene_IDs (common name) for selected FILE/all files in selected FOLDER'
	print

else:
	if os.path.isfile(sys.argv[1]):
		targetf= [sys.argv[1]]
		parent= ''
	else:
		targetf=[]
		for f in os.listdir(sys.argv[1]):
			if f[-3:]=='.fa':
				targetf.append(f)
		parent= sys.argv[1]
		if parent[-1]!='/':
			parent+='/'				

	for t in targetf:
		outfile= t.split('.')[0]+'.map'
		grep= sub.Popen('grep ">" '+parent+t, shell=True, stdout=1)
		lines= grep.communicate()[0].split('\n')
		f=open(outfile,'w')
		for line in lines:
			if line!='' and len(line.split())>=2:
				x=line.split()
				NAME= x[0][1:]
				PAC= x[1].split('PACid:')[1]
				f.write('PAC:'+PAC+'\t'+NAME+'\n')
		f.close()
				
			 
