#!/usr/bin/python

import sys
import subprocess as sub
sub.PIPE=1

Ncores= sub.Popen('grep processor /proc/cpuinfo | wc -l', shell=True, stdout=1)
Nthreads= int(Ncores.communicate()[0])-1 #one core used to run current process

EXpath= 'Exonerate' #'/home/alessia/programs/exonerate-2.2.0-x86_64/bin/exonerate'
N= 10
ryo= '"\ntarget \tSTART %tab END %tae LENGTH %tal\n\n"'

if len(sys.argv)!=4:
	print '\nusage: <query.fa> <target.fa> <alignment model (#)>'
	print '\ndefault N hits:\t\t'+str(N)
	print '\nmodel #:'
	print '\t1. protein2dna:bestfit\tN hits up to sequences in query - sometimes more due to refinement process'
	print '\t\t\t\t(NOTE: full query coverage, i.e. all alignments to a query have same length)'
	print '\t2. protein2dna\t\tdefault N hits'
	print '\t3. protein2genome\tdefault N hits'
	print
else:
	query= sys.argv[1]
	target= sys.argv[2]
	m= int(sys.argv[3])

	if m not in [1,2,3]:
		print 'invalid model # provided: '+str(m)
	else:
		if m==1:
			model='protein2dna:bestfit'
			E= True
			refine= 'none'
			x= sub.Popen('grep ">" '+target+' | wc -l', shell=True, stdout=1)
			N= int(x.communicate()[0])
		elif m==2:
	 		model='protein2dna'
			E=False
			refine= 'full' #aligned sequences both fully refined (possibly extended)
		else: #m==3
 			model='protein2genome'
			refine= 'region --refineboundary 100' #aligned sequences refined only in alignment region + boundary
			E=False
		Q= 'protein'
		T= 'dna'
	
		command= ''+EXpath+' -q '+query+' -Q '+Q+' -t '+target+' -T '+T+' -E '+str(E)+' -m '+model+' -n '+str(N)+' -s 1 --refine '+refine+' --alignmentwidth 100 --ryo '+ryo+''  
		#x=sub.Popen('echo '+command, shell=True, stdout=1) #check how the command passed to the shell looks like
		x=sub.Popen(command, shell=True, stdout=1)
		OUT=x.communicate()[0]
	
		print OUT #prints to standard output
		Nout= len(OUT.split('vulgar'))-1
		if Nout>0:
			print '\ntotal hits count:\t'+str(Nout)+'\n'

