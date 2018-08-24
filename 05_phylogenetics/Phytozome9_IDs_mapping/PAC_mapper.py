#!/usr/bin/python

import sys
import os
import subprocess as sub
sub.PIPE=1

default= '/linuxhome/alessia2/PROJECT02_bZIPs_tree/Phytozome9_protein_original_IDs/' #note: final '/'!

if len(sys.argv)!=2:
	print 'usage: <file OR folder> (write "default" or substring for default folder - see below)'	
	print
	print 'maps PAC_IDs to gene_IDs (commond name) for selected FILE/all files in selected FOLDER'
	print
	print 'default folder:'
	print default
	print
else:
	if len("default".split(sys.argv[1].lower()))==1:
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
	else:
		targetf= []
		for f in os.listdir(default):
			if f[-3:]=='.fa':
				targetf.append(f)
		parent= default

	for t in targetf:
		outfile= t.split('_')[0]+'.map'
		grep= sub.Popen('grep ">" '+parent+t, shell=True, stdout=1)
		lines= grep.communicate()[0].split('\n')
		f=open(outfile,'w')
		for line in lines:
			if line!='':
				x=line.split('|PACid:')
				NAME= x[0].split()[0][1:]
				PAC= x[1].split()[0]
				f.write('PAC:'+PAC+'\t'+NAME+'\n')
		f.close()
				
			 
