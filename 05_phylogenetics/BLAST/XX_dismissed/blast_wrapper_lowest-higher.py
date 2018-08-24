#!/usr/bin/python

import sys
import subprocess as sub
sub.PIPE=1

#Ncores= sub.Popen('grep processor /proc/cpuinfo | wc -l', shell=True, stdout=1)
Nthreads= 1 #int(Ncores.communicate()[0])-1 #one core used to run current process

BLAST= 'blastp' #'tblastn'
MIN= True	#if true, finds minimum score for each hit among all pairwise comparisons (might be useful for closely related seqs, ex. bZIPs)
E= 1e-10	#max. e-value 

#			S	
#Ath vs Ath (full):	719	
#Osa vs Osa (full): 	721	
	

if len(sys.argv)!=4:
	print
	print '###########################################################'
	print '# BLASTP HITS with score > SCORE (for e-value check code) #'
	print '# >usage: blastP <subject> <query> <SCORE>   	         #'
	print '# also available MIN option (check code)                  #'
	print '# sometimes useful for closely related/ all vs all seqs   #'
	print '# (ex. Ath bZIPs vs Ath bZIPs)			         #'
	print '###########################################################'
	print '\npre-comparison filtering:'
	print '- maximum e-value: '+str(E)
	print '- BLAST: '+BLAST
	if MIN:
		print '- MIN option ON!'
	print
else:
	subj= sys.argv[1]
	q= sys.argv[2]
	S= int(sys.argv[3])
	out= '.'.join(subj.split('/')[-1].split('.')[:-1])+'__'+'.'.join(q.split('/')[-1].split('.')[:-1])+'.out'
	if BLAST=='blastp':
		BLAST2= BLAST
	elif BLAST=='tblastn':
		BLAST2= 'tblastx'
	elif BLAST=='tblastx':
		BLAST2=='tblastn'

	print
	print '#GENERAL SETTINGS'
	print '\tBLAST:',BLAST
	print '#pre-comparison filtering:'
	print '\tmin score:',S
	print '\tmax e-value:',E
	if MIN:
		print '\t***MIN option ON!***'
	#BLAST 

	print 'starting blast subj. vs query..'

	blast= sub.Popen(BLAST+' -subject '+subj+' -query '+q+' -evalue '+str(E)+' -outfmt "6 sseqid score evalue" -num_threads '+str(Nthreads),shell=True, stdout=1)
	hits= blast.communicate()[0].split('\n')

	print 'DONE'

	#filter hits (score and e-value)

	HITS= {}
	for h in hits:
		x= h.split('\t')
		if len(x)!=3:	
			continue	
		ID=x[0]
		score= int(x[1])
		if x[2][0]=='e':
			evalue= float('1'+x[2])
		else:
			evalue= float(x[2])
		if score >= S: #e-value filtering already in blast, not necessary here
			if HITS.get(ID,0)==0 or (MIN and (HITS[ID][1]<evalue or (HITS[ID][1]==evalue and HITS[ID][0]>score))) or (not MIN and (HITS[ID][1]>evalue or (HITS[ID][1]==evalue and HITS[ID][0]<score))): 
				HITS[ID]=(score,evalue) 
	hits= None #free memory

	if len(HITS.keys())==0:
		print 'NO SIGNIFICANT HITs!'

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
