#!/usr/bin/python

import sys
import time
import subprocess as sub
sub.PIPE=1

def count(ITER):
	n=0
	for a in ITER:
		if a:#True
			n+=1
	return n

def all(ITER):
	A=True
	for a in ITER:
		if not a:
			A=False
			break
	return A

#NOTE: multithread version will allow more hits for same N compared to single core!
#only SCORE ensures same number of hits between single/multiprocess run

Ncores= int(sub.Popen('grep processor /proc/cpuinfo | wc -l', shell=True, stdout=1).communicate()[0]) #one core used to run main process
EXpath= 'Exonerate' #'/home/alessia/programs/exonerate-2.2.0-x86_64/bin/exonerate'
#ryo= '"\ntarget \tSTART %tab END %tae LENGTH %tal\n\n"' #reports custom information in more readable format

if len(sys.argv)!=4:
	print '\nusage: <query.fa> <target.fa> <alignment model (1|2)>'
	print '\ncompares each individual query to each individual target (1 hit allowed) and outputs results sorted by score'
	print '(at most N query x N target sequences)'
	print '\nmodel #:'
	print '\t1. protein2dna:bestfit\t(tries to achieve full query coverage, i.e. hits more likely to include start/stop codon)'
	print '\t2. protein2dna\t\t(standard, local alignment - usually finds more hits than previous, but shorter)'
	print
else:
	query= sys.argv[1]
	target= sys.argv[2]
	m= int(sys.argv[3])
	Nquery= int(sub.Popen('grep ">" '+query+' | wc -l', shell=True, stdout=1).communicate()[0])
	if m not in [1,2]:
		print 'invalid model # provided: '+str(m)
	else:	#MODEL CHECK
		if m==1:
			model='protein2dna:bestfit'
			E= True
			refine= 'none'

		else: #m==2:
	 		model='protein2dna'
			E=False
			refine= 'full' #aligned sequences both fully refined (possibly extended)
		Q= 'protein'
		T= 'dna'

		outfile= '.'.join(target.split('.')[:-1])+'.vs.'+'.'.join(query.split('/')[-1].split('.')[:-1])+'.ex'+str(m)
		
		#get N targets and their position in target file
		IDlines= sub.Popen('grep -n ">" '+target+'',shell=True,stdout=1).communicate()[0].split('\n')
		nlines= sub.Popen('cat '+target+' | wc -l',shell=True,stdout=1).communicate()[0]
		nLIST=[]
		IDcheck=[] #store target names for final check
		for line in IDlines:
			if len(line.split(':'))>1:
				nLIST.append(int(line.split(':')[0])) #lines containing sequence IDs in the target file
				IDcheck.append(':'.join(line.split(':')[1:])[1:-1]) #excludes '>' and '\n'
		nLIST.append(int(nlines)+1) #index of last line
		nID=len(nLIST)-1
		#MULTITHREAD PROCESSING of query against individual target sequences
		JOBS={}
		for ID in range(nID):
			while count(map(lambda x: JOBS[x].poll()==None,JOBS.keys()))==Ncores-1: #max threads allowed at one time
				#print 'now sleeping1'
				time.sleep(1)
			print '\nfinding hits for target',ID+1,'of',nID,'.. (',IDcheck[ID],')'
			tmp= 'tmp'+str(ID)+'.tmp'
			sub.call('sed -n "'+str(nLIST[ID])+','+str(nLIST[ID+1]-1)+' p" '+target+' > '+tmp+'',shell=True)
			JOBS[ID]=sub.Popen(''+EXpath+' -q '+query+' -Q '+Q+' -t '+tmp+' -T '+T+' -E '+str(E)+' -m '+model+' -n 1 -s 1 --refine '+refine+' --alignmentwidth 200 && rm '+tmp+'',shell=True,stdout=1)
		while not all(map(lambda x: JOBS[x].poll()!= None,JOBS.keys())):#waits for all the processes to be concluded
			#print 'now sleeping2'
			time.sleep(1)
		#collect output when all processes are concluded
		OUT=''
		for ID in JOBS.keys():
			OUT+= JOBS[ID].communicate()[0]

		#OUTPUT PROCESSING
		#check missing IDs
		BAD=[]
		for ID in IDcheck:
			if OUT.find(ID)==-1: #failure to find ID
				BAD.append(ID)
		#results filtering
		split1='C4 Alignment:'
		split2='-- completed exonerate analysis'
		OUT= OUT.split(split1)
		if len(OUT)>1:
			OUT=OUT[1:]
			SDB=[]
			for i in range(len(OUT)):
				q=OUT[i]
				score=int(q.split('Raw score: ')[1].split()[0])
				SDB= SDB+[(score,q.split(split2)[0])]
			SDB.sort(reverse=True)
			OUT=''
			for q in SDB:
				OUT+= split1+q[1]
			f=open(outfile,'w')
			f.write(OUT)
			f.write('\n*************************************************\n')
			f.write('#expected max N:\t'+str(Nquery*nID)+' ('+str(Nquery)+' queries x '+str(nID)+' targets)\n')
			f.write('#retrieved hits:\t'+str(len(SDB))+'\n')
			#print OUT
			print
			print '*************************************************'
			print '\nresults stored in:\t'+outfile
			print '\nexpected max N:\t'+str(Nquery*nID)+' ('+str(Nquery)+' queries x '+str(nID)+' targets)'
			print 'retrieved hits:\t'+str(len(SDB))
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
			
