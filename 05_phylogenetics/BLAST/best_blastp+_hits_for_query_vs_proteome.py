#!/usr/bin/python

import sys
import os
import subprocess as sub
sub.PIPE=1

#Ncores= sub.Popen('grep processor /proc/cpuinfo | wc -l', shell=True, stdout=1)
Nthreads= 1 #int(Ncores.communicate()[0])-1 #one core used to run current process

S= 50		#min. score (query1 initial filtering)
E= 1e-5		#max. e-value (query1 initial filtering)
N= 2

#			S	E	N
#Ath+Osa vs algae:	30	1e-5	5 (2 for stringent)
#Ath+Osa vs ancient: 	50	1e-5	2
#Ath SbZIPs vs Osa:	200	1e-10	2

if len(sys.argv)!=4:
	print
	print '######################################################################################'
	print '# finds QUERY BLASTP HITS and re-blasts against proteome from which query comes from #'
	print '# to keep only hits that have at least 1 query sequence within first N results       #'
	print '# >usage: blastP <subject> <query> <proteome (incl. query)> (fasta format)           #'
	print '# (e.g. <algae proteins> <C+S bZIPs Ath+Osa> <proteome Ath+Osa (incl. C+S bZIPs)>    #'
	print '######################################################################################'
	print
	print 'pre-filtering:'
	print '- minimum score: '+str(S)
	print '- maximum e-value: '+str(E)
	print
	print 'hits included in final output if:'
	print '- reverse Blast includes a query sequence in the first '+str(N)+' hits'
	print
else:
	subj= sys.argv[1]
	q= sys.argv[2]
	p= sys.argv[3]
	out= '.'.join(subj.split('/')[-1].split('.')[:-1])+'__'+'.'.join(q.split('/')[-1].split('.')[:-1])+'__'+'.'.join(q.split('/')[-1].split('.')[:-1])+'.rev'+str(N)+'.out'
	tmp= '.'.join(subj.split('/')[-1].split('.')[:-1])+'__'+'.'.join(q.split('/')[-1].split('.')[:-1])+'__'+'.'.join(q.split('/')[-1].split('.')[:-1])+'.rev'+str(N)+'.tmp'
	
	print
	print '#GENERAL SETTINGS'
	print '#pre-comparison filtering:'
	print '\tmin score:',S
	print '\tmax e-value:',E
	print '#hits kept if:'
	print '\trev. Blast contains query in first',N,'hits'
	print

	#BLAST QUERY

	print 'starting blast subj. vs query..'

	blastq= sub.Popen('blastp -subject '+subj+' -query '+q+' -evalue '+str(E)+' -outfmt "6 sseqid score evalue" -num_threads '+str(Nthreads),shell=True, stdout=1)
	hitsq= blastq.communicate()[0].split('\n')

	print 'DONE'

	#filter hits (score and e-value)

	HITS= {}
	for h in hitsq:
		x= h.split('\t')
		if len(x)!=3:	
			continue	
		ID=x[0]
		score= int(x[1])
		if x[2][0]=='e':
			evalue= float('1'+x[2])
		else:
			evalue= float(x[2])
		if score >= S:
			if HITS.get(ID,0)==0:
				HITS[ID]=(score,evalue)
			elif HITS[ID][1]>evalue or (HITS[ID][1]==evalue and HITS[ID][0]<score):
				HITS[ID]=(score,evalue) 
	hitsq= None #free memory

	if len(HITS.keys())==0:
		print 'NO SIGNIFICANT HITs!'
	else:
		#extract hits sequences for inverse blast
		f=open(tmp,'w')
		APPEND=False
		for line in open(subj).readlines():
			if len(line)>1:
				if line[0]=='>':
					if HITS.get(line[1:].split()[0],0)!=0:
						f.write(line)
						APPEND=True
					else:
						APPEND=False
				elif APPEND:
					f.write(line)
		f.close()

		#REVERSE-BLAST PROTEOME

		print 'starting reverse-blast hits. vs proteome..'

		blastp= sub.Popen('blastp -subject '+p+' -query '+tmp+' -evalue '+str(E)+' -max_target_seqs '+str(N)+' -outfmt "6 qseqid sseqid" -num_threads '+str(Nthreads),shell=True, stdout=1)
		hitsp= blastp.communicate()[0].split('\n')

		print 'DONE'
		
		os.remove(tmp)

		#get reverse hits
	
		REV={}
		for h in hitsp:
			x= h.split('\t')	
			if len(x)!=2:
				continue
			REV[x[0].split()[0]]=REV.get(x[0].split()[0],[])+[x[1].split()[0]]
		hitsp= None #free memory

		#get query IDs

		Qseqs={}
		for line in open(q).readlines():
			if len(line)>1 and line[0]=='>':
				Qseqs[line[1:].split()[0]]=1 #removes initial '>' and everything after 1st space 
								#(apparently blast does the same to sequence names)

		#filter hits (remove hits with no query IDs among reverse hits)
	
		for r in REV.keys():
			OK=False
			for i in REV[r]:
				if Qseqs.get(i,0)!=0:
					OK=True
					break
			if not OK:
				HITS.pop(r)

		#output filtered blast out

		if len(HITS.keys())==0:
			print 'NO REVERSE HITs!'
		else:

			H= HITS.keys()
			H.sort()
			f=open(out,'w')
			for h in H:
				if len(h.split('|'))==1:
					ID= h.split()[0]
				else:
					ID= '|'.join(h.split('|')[:-1])
				values = '\t%d\t%.2e' % (HITS[h][0], HITS[h][1])
				f.write(ID+values+'\n')
			f.close() 

			print '\nfiltered hits are in',out
