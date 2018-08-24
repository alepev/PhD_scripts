#!/usr/bin/python

import sys
import subprocess as sub
sub.PIPE=1

#Ncores= sub.Popen('grep processor /proc/cpuinfo | wc -l', shell=True, stdout=1)
Nthreads= 1 #int(Ncores.communicate()[0])-1 #one core used to run current process

BLAST= 'tblastn'

S= 50		#min. score (query1 initial filtering)
E= 1e-5		#max. e-value (query1 initial filtering)
Eratio=10	#hit query1 must have e-value no larger than hit query2 e-value * Eratio
Sgap=20		#hit query1 must have score no smaller than hit query2 score - Sgap

#			S	E	Eratio	Sgap
#Ath vs Brassicales: 	100	1e-10	10	10
#Ath vs dicots:		90	1e-5	10	20
#Osa vs monocots:	"	"	"	"	
#Ath/Osa vs ancient:	50	1e-5	10	20

if len(sys.argv)!=4:
	print
	print '#####################################################################'
	print '# finds QUERY1 BLASTP HITS with better e-value or score than query2 #'
	print '# >usage: blastP <subject> <query1> <query2> (fasta format)         #'
	print '# e.g. <dicot proteins> <C bZIPs from Ath> <non-C bZIPs from Ath>   #'
	print '#####################################################################'
	print '\npre-comparison filtering:'
	print '- minimum score: '+str(S)
	print '- maximum e-value: '+str(E)
	print
	print 'borderline cases (included in final output):'
	print '- e-value1 <= '+str(Eratio)+' * e-value2' 
	print '- equal e-values with score1 >= (score2 - '+str(Sgap)+')'
	print
else:
	subj= sys.argv[1]
	q1= sys.argv[2]
	q2= sys.argv[3]
	out= '.'.join(subj.split('/')[-1].split('.')[:-1])+'__'+'.'.join(q1.split('/')[-1].split('.')[:-1])+'__'+'.'.join(q2.split('/')[-1].split('.')[:-1])+'.out'
	if BLAST=='blastp':
		BLAST2= BLAST
	elif BLAST=='tblastn':
		BLAST2= 'tblastx'
	elif BLAST=='tblastx':
		BLAST2=='tblastn'

	print
	print '#GENERAL SETTINGS'
	print '#pre-comparison filtering:'
	print '\tmin score:',S
	print '\tmax e-value:',E
	print '#borderline cases kept if:'
	print '\te-value1 <= e-value2 *',Eratio
	print ' OR (if same e-value)'
	print '\tscore1 >= ( score2 -',Sgap,')'
	print

	#BLAST QUERY1 

	print 'starting blast subj. vs query1..'

	blast1= sub.Popen(BLAST+' -subject '+subj+' -query '+q1+' -evalue '+str(E)+' -outfmt "6 sseqid score evalue" -num_threads '+str(Nthreads),shell=True, stdout=1)
	hits1= blast1.communicate()[0].split('\n')

	print 'DONE'

	#filter hits (score and e-value)

	HITS= {}
	for h in hits1:
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
			if HITS.get(ID,0)==0:
				HITS[ID]=(score,evalue)
			elif HITS[ID][1]>evalue or (HITS[ID][1]==evalue and HITS[ID][0]<score):
				HITS[ID]=(score,evalue) 
	hits1= None #free memory

	if len(HITS.keys())==0:
		print 'NO SIGNIFICANT HITs!'

	else:

		#BLAST QUERY2

		print 'starting blast subj. vs query2..'
	
		blast2= sub.Popen(BLAST2+' -subject '+subj+' -query '+q2+' -evalue '+str(E)+' -outfmt "6 sseqid score evalue" -num_threads '+str(Nthreads),shell=True, stdout=1)
		hits2= blast2.communicate()[0].split('\n')

		print 'DONE'
	
		#filter hits (comparison with hits query2)
	
		for h in hits2:
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
				if HITS.get(ID,0)!=0 and ((HITS[ID][1]> evalue*Eratio) or (HITS[ID][1]==evalue and HITS[ID][0]<score-Sgap)):				
					HITS.pop(ID)
		hits2= None #free memory

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
