#!/usr/bin/python

import sys
import time
import subprocess as sub
sub.PIPE=1

def all(ITER):
	A=True
	for a in ITER:
		if not a:
			A=False
			break
	return A

#NOTE: multithread version will allow more hits for same N compared to single core!
#only SCORE ensures same number of hits between single/multiprocess run

#PROBLEM: detected malfunctioning of protein2dna:bestfit with multithread version!

Ncores= sub.Popen('grep processor /proc/cpuinfo | wc -l', shell=True, stdout=1)
Ncores= 1 #int(Ncores.communicate()[0])-1 #one core used to run current process

EXpath= 'Exonerate' #'/home/alessia/programs/exonerate-2.2.0-x86_64/bin/exonerate'
#ryo= '"\ntarget \tSTART %tab END %tae LENGTH %tal\n\n"' #reports custom information in more readable format

if len(sys.argv)<4 or len(sys.argv)>5:
	print '\nusage: <query.fa> <target.fa> <alignment model (#)> <N (optional)>'
	print '\ndefault N hits (if not provided by user): N QUERY SEQUENCES (assumes 1 match per query)'
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
	if len(sys.argv)==5:
		N=int(sys.argv[4])
		Nstring= '(user-provided)'
	else:
		x= sub.Popen('grep ">" '+target+' | wc -l', shell=True, stdout=1)
		N= int(x.communicate()[0])
		Nstring= '(N target sequences)'

	if m not in [1,2,3]:
		print 'invalid model # provided: '+str(m)
	else:	#MODEL CHECK
		if m==1:
			model='protein2dna:bestfit'
			E= True
			refine= 'none'

		elif m==2:
	 		model='protein2dna'
			E=False
			refine= 'full' #aligned sequences both fully refined (possibly extended)
		else: #m==3
 			model='protein2genome'
			E=False
			refine= 'region --refineboundary 100' #aligned sequences refined only in alignment region + boundary
		Q= 'protein'
		T= 'dna'

		outfile= '.'.join(target.split('.')[:-1])+'.vs.'+'.'.join(query.split('/')[-1].split('.')[:-1])+'.ex'

		#load target IDs to check later if all of them got a hit
		IDcheck={}
		for line in open(target).readlines():
			if line!='' and line[0]=='>':
				IDcheck[line[1:-1]]=1 #removes '>' and '\n'
		#COMMAND EXECUTION
		command= ''+EXpath+' -q '+query+' -Q '+Q+' -t '+target+' -T '+T+' -E '+str(E)+' -m '+model+' -n '+str(N)+' -s 1 --refine '+refine+' --alignmentwidth 100'
		#add option for customized output information: --ryo '+ryo+''  
		if Ncores<2:
			JOB= sub.Popen(command, shell=True, stdout=1)
			OUT=JOB.communicate()[0]
		else:
			JOBS={}
			for i in range(Ncores):
				#command debug - check how the command passed to the shell looks like
				#JOBS[i]= sub.Popen('echo "'+command+' --targetchunkid '+str(i+1)+' --targetchunktotal '+str(Ncores)+'"', shell=True, stdout=1) 
				JOBS[i]= sub.Popen(command+' --targetchunkid '+str(i+1)+' --targetchunktotal '+str(Ncores)+'', shell=True, stdout=1)
			while not all(map(lambda x: JOBS[x].poll()==0,JOBS.keys())): #WAITS FOR ALL THE PROCESSES TO BE CONCLUDED
				time.sleep(1) #checking time (in seconds)
			OUT=''
			Nout=0
			for i in range(Ncores):
				OUT+=JOBS[i].communicate()[0]
		#OUTPUT
		#files + command info
		split1='C4 Alignment:'
		split2='-- completed exonerate analysis'
		OUT= OUT.split(split1)
		#results filtering
		if len(OUT)>1:
			f=open(outfile,'w')
			SDB=[]
			for q in OUT[1:]:
				score=int(q.split('Raw score: ')[1].split()[0])
				SDB= SDB+[(score,q.split(split2)[0])]
			SDB.sort(reverse=True)
			OUT=''
			if len(SDB)>N:
				SDB=SDB[:N] #output only first N results
			BAD=[]
			for ID in IDcheck.keys():
				if SDB.find(ID)==-1: #failure to find ID
					BAD.append(ID)
			for q in SDB:
				OUT+= split1+q[1]
			f.write(OUT)
			f.write('\n*************************************************\n')
			f.write('#SETTINGS: max.N == '+str(N)+' '+Nstring+'\n')
			#print OUT
			print
			print '*************************************************'
			print '\nresults stored in:\t'+outfile
			print '\nSETTINGS:\t max.N == '+str(N)+' '+Nstring
			if len(SDB)<N:
				f.write('\n*** WARNING: N hits ('+str(len(SDB))+') less than max. N! ***\n')
				print '\n*** WARNING: N hits ('+str(len(SDB))+') less than max. N! ***'
			if len(BAD)>0:
				f.write('\n*** WARNING: the following target sequence(s) didn\'t produce any result! ***\n')
				print '\n*** WARNING: the following target sequence(s) didn\'t produce any result! ***'
				BAD.sort()
				for ID in BAD:
					f.write(ID+'\n')
					print ID
			f.close()
		else:
			print '\n***no result found!***'
		print
			
