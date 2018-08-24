#!/usr/bin/python

import sys
import copy

if len(sys.argv) not in [2,3]:
	print
	print 'usage: <BLAST output> <SCORE threshold> (OPTIONAL: default ALL hits)'
	print '(NOTE: assumes STANDARD blast legacy / blast+ output from search with -DB OPTION)'
	print
	print '	OUTPUT:'
	print '	 - hits sequence limits file (MERGED overlapping + original hits commented)'
	print '	 - hits aligned regions fasta file (original hits only)'
	print '	(MERGED aligned region is calculated from overlapping hits; e-value/score from best hit)'
	print 

#ABOUT BLAST INDEXING:
#index system starting with 1, end inclusive  
#ex. sequence 1:4:
#     A C G T
#     1 2 3 4
#this notation is mantained also in the output files.
#to obtain corresponding string: seq[start-1:end] 
#(inclusive/exclusive, respectively)
#ex. the previous becomes sequence[0:4]

else:
	blast= sys.argv[1]
	if len(sys.argv)==3:
		score= float(sys.argv[2])
		ALL=False
	else:
		ALL=True

	fastaout= '.'.join(blast.split('.')[:-1])+'.out.fa'#'.'.join(blast.split('/')[-1].split('.')[:-1])+'.out.fa'
	limitout= '.'.join(blast.split('.')[:-1])+'.out.limits.txt'#'.'.join(blast.split('/')[-1].split('.')[:-1])+'.out.limits.txt'
	
	#extract blast results
	HITS= {}
	x= '\n'.join(open(blast).readlines()).split('Query= ')
	for query in x[1:]:
		QUERY= query.split('\n')[0].rstrip()
		y= query.split('> ') #-db output option
		if len(y)==1:
			y= query.split('Subject= ') #-subject output option
			if len(y)==1:
				print '*** no hits found for',QUERY,'***'
				continue #NO HIT FOR THAT QUERY
		for subj in y[1:]: #many hits can be reported for same subject when long target sequences 
			HIT= subj.split('\n')[0].rstrip()
			y2= subj.split('Score = ')
			for hit in y2[1:]:
				S= float(hit.split(' bits')[0])
				E= float(hit.split('Expect')[1].split(',')[0].split()[-1]) #normally "Expect =", however when multiple hits in -db option -> "Expect(N) = "
				if not ALL and S<score:
					continue #HIT SCORE BELOW THRESHOLD
				z=hit.split('Sbjct')
				START= int(z[1].split('\n')[0].lstrip().split()[0])
				if len(z)>2:
					STOP= int(z[-1].split('\n')[0].rstrip().split()[-1])
					SEQ=''
					for seq in z[1:]:
						SEQ+= seq.split('\n')[0].lstrip().rstrip().split()[1]
				else:
					STOP= int(z[1].split('\n')[0].rstrip().split()[-1])
					SEQ= z[1].split('\n')[0].lstrip().rstrip().split()[1]
				if START<STOP:
					STRAND='+'
				else:
					START,STOP=STOP,START
					STRAND='-'
				SEQ=SEQ.replace('-','').upper()
				HITS[HIT]= HITS.get(HIT,[])+[(E,S,START,STOP,STRAND,SEQ,QUERY)]

	#verify overlap between stored values & merge 
	if len(HITS.keys())>0:
		Mall={}
		for HIT in HITS.keys():
			M1,M2=[],[]
			for hit in HITS[HIT]:
				if hit[4]=='+': #STRAND
					M1+=[(hit[2],hit[3],hit[0],hit[1])] #START,STOP,E,S
				else:
					M2+=[(hit[2],hit[3],hit[0],hit[1])] #START,STOP,E,S
			for LIST in [M1,M2]:
				if LIST!=[]:
					if LIST==M1:
						STRAND='+'
					else:
						STRAND='-'
					MERGED=copy.copy(LIST)
					MERGED.sort()
					CHECK=True
					done=0
					while len(MERGED)>1 and CHECK:
						CHECK=False
						M=MERGED[:done]
						for i in range(done,len(MERGED)-1):
							A,B=MERGED[i][0],MERGED[i][1]
							E1,S1=MERGED[i][2],MERGED[i][3]
							C,D=MERGED[i+1][0],MERGED[i+1][1]
							E2,S2=MERGED[i+1][2],MERGED[i+1][3]
							#print A,B,C,D
							if B<C:
								M+=[(A,B,E1,S1)]
								if i+1==len(MERGED)-1: #last iteration so no further comparisons
									M+=[(C,D,E2,S2)]
								#print 'B<C'
							else: #overlap! E,S to save best score among merged hits
								if E1<E2:
									E,S=E1,S1
								elif E1>E2:
									E,S=E2,S2
								else:#E1==E2
									E=E1
									if S1>=S2:
										S=S1
									else:
										S=S2	
								if B<=D:	
									M+=[(A,D,E,S)]	
									#print 'C<B<D'
								else:
									M+=[(A,B,E,S)]
									#print 'B>D'
								if i+1<len(MERGED)-1:
									M+=MERGED[i+2:]
									CHECK=True
									done=i
								break #restart iteration with updated MERGED
						MERGED=copy.copy(M)
					for m in MERGED:
						#print m
						Mall[HIT]=Mall.get(HIT,[])+[(m[2],-m[3],m[0],m[1],STRAND)] #E,-S,START,STOP + STRAND (-S to allow sorting in decreasing order as for E)
			m= Mall[HIT]
			m.sort()
			Mall[HIT]=m

		#output summary
		f= open(fastaout,'w')
		l= open(limitout,'w')
		l.write('#ID\tSTART\tSTOP\tSTRAND\tSEQLEN\tE-VALUE\tSCORE\tQUERY\n')
		print
		print '#ID\tSTART\tSTOP\tSTRAND\tSEQLEN\tE-VALUE\tSCORE\tQUERY'
		SORTED=[]
		for HIT in HITS.keys():
			scores=HITS[HIT]
			scores.sort()
			HITS[HIT]=scores
			SORTED+= [(scores[0][0],-scores[0][1],HIT)] #best E,S for hit (-S to allow sorting in decreasing order as for E)
		SORTED.sort()
		for sorted in SORTED:
			HIT=sorted[2]
			for m in Mall[HIT]:
				E,S,START,STOP,STRAND=m[0],-m[1],m[2],m[3],m[4]
				l.write(HIT+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t'+str(E)+'\t'+str(S)+'\tMERGED\n')
				print HIT+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t'+str(E)+'\t'+str(S)+'\tMERGED'
			for score in HITS[HIT]:
				E=score[0]
				S=score[1]
				START= score[2]
				STOP= score[3]
				STRAND= score[4]
				SEQ= score[5]
				QUERY= score[6]
				f.write('>'+HIT+' | '+str(START)+':'+str(STOP)+' ['+STRAND+'] ('+str(STOP-START+1)+'nt)- blast vs '+QUERY+' (E'+str(E)+',S'+str(S)+')\n')
				f.write(SEQ+'\n')
				#output preceded by '#' to allow further processing of MERGED hits only (but original remain visible for inspection)
				l.write('#'+HIT+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t'+str(E)+'\t'+str(S)+'\t'+QUERY+'\n')
				print '#'+HIT+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t'+str(E)+'\t'+str(S)+'\t'+QUERY
				#print SEQ	
		f.close()
		l.close()

		print
		print '*********************************************************'
		print 'check output files for results:'
		print ' - aligned sequences limits:',limitout
		print ' - aligned sequences fasta:', fastaout
		print
		print 'NOTE! SEQLEN refers to DNA SEQUENCE in tblastn & tblastx (even if reported alignment is protein)'
		print

	else:
		print
		print '*** no blast hit retrieved! ***'
		print
	
