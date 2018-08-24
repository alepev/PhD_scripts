#!/usr/bin/python

import sys,copy,os
import subprocess as sub

if len(sys.argv)!=2:
	print 'usage: runs WEEDER with target files from folder <target_folder>'
	print 'wrapper <target_folder>'
else:
	script= '/home/alessia/programs/Weeder1.4.2/weederlauncher.out'
	target_folder= sys.argv[1]
	if target_folder[-1]!='/':
		target_folder+='/'
	target_list=os.listdir(target_folder)
	for target in target_list:
		if target[-4:]=='.txt':
#			print target_folder+target
			sub.call('echo '+script+' '+target_folder+target+' AT medium S T20',shell=True)
			sub.call(script+' '+target_folder+target+' AT medium S T20',shell=True)
			sub.call('echo "#======================================================================"',shell=True)
