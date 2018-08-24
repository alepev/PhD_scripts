#!/usr/bin/python

import sys

if len(sys.argv) not in [3,4] or not sys.argv[2].isdigit() or (len(sys.argv)==4 and sys.argv[3]!='1'):
	print
	print 'usage: <BLAST output> <e-value threshold (0 for all)> <1 for ALL SEQS (optional)>'
	print '(NOTE: assumes standard unformatted blast legacy / blast+ output from search with -db option)'
	print
	print '	OUTPUT:'
	print '	 - hits sequence limits file'
	print '	 - hits aligned regions fasta file'
	print '	NOTE: for each hit only aligned region with best e-value is kept'
	print '	(to keep ALL alignments use option 1)'
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
	evalue= sys.argv[2]

	if evalue=='0':
		ALL=True
	else:
		evalue=float(evalue)
		ALL=False

	if len(sys.argv)==4 and sys.argv[3]=='1':
		MULTI=True
		fastaout= '.'.join(blast.split('/')[-1].split('.')[:-1])+'_ALL.out.fa'
		limitout= '.'.join(blast.split('/')[-1].split('.')[:-1])+'_ALL.out.limits.txt'
	else:
		MULTI=False
		fastaout= '.'.join(blast.split('/')[-1].split('.')[:-1])+'_BEST.out.fa'
		limitout= '.'.join(blast.split('/')[-1].split('.')[:-1])+'_BEST.out.limits.txt'

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
				if E>evalue and not ALL:
					continue #HIT E-VALUE ABOVE THRESHOLD
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
				if MULTI:
					HITS[HIT]= HITS.get(HIT,[])+[(E,S,START,STOP,STRAND,SEQ,QUERY)]
				else:
					if HITS.get(HIT,0)==0 or HITS[HIT][0]>E or (HITS[HIT][0]==E and HITS[HIT][1]<S):
						HITS[HIT]= [E,S,START,STOP,STRAND,SEQ,QUERY] 

	#output summary
	
	if len(HITS.keys())>0:
		f= open(fastaout,'w')
		l= open(limitout,'w')
		l.write('#ID\tSTART\tSTOP\tSTRAND\tSEQLEN\tE-VALUE\tSCORE\tQUERY\n')
		print
		print '#ID\tSTART\tSTOP\tSTRAND\tSEQLEN\tE-VALUE\tSCORE\tQUERY'
		if MULTI:
			SORTED=[]
			for hit in HITS.keys():
				scores=HITS[hit]
				scores.sort()
				HITS[hit]=scores
				SORTED+= [(scores[0][0],-scores[0][1],hit)] #best E,S for hit (-S to allow sorting in decreasing order as for E)
			SORTED.sort()
			for sorted in SORTED:
				hit=sorted[2]
				for score in HITS[hit]:
					E=score[0]
					S=score[1]
					START= score[2]
					STOP= score[3]
					STRAND= score[4]
					SEQ= score[5]
					QUERY= score[6]
					f.write('>'+hit+' | '+str(START)+':'+str(STOP)+' ['+STRAND+'] ('+str(STOP-START+1)+'nt)- blast vs '+QUERY+' (E'+str(E)+',S'+str(S)+')\n')
					f.write(SEQ+'\n')
					l.write(hit+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t'+str(E)+'\t'+str(S)+'\t'+QUERY+'\n')
					print hit+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t'+str(E)+'\t'+str(S)+'\t'+QUERY
					#print SEQ
		else:
			SORTED=[]
			for hit in HITS.keys():
				SORTED+= [(HITS[hit][0],-HITS[hit][1],hit)] #best E,S for hit (-S to allow sorting in decreasing order as for E)
			SORTED.sort()
			for sorted in SORTED:
				hit=sorted[2]
				score= HITS[hit]
				E=score[0]
				S=score[1]
				START= score[2]
				STOP= score[3]
				STRAND= score[4]
				SEQ= score[5]
				QUERY= score[6]
				f.write('>'+hit+' | '+str(START)+':'+str(STOP)+' ['+STRAND+'] ('+str(STOP-START+1)+'nt)- blast vs '+QUERY+' (E'+str(E)+',S'+str(S)+')\n')
				f.write(SEQ+'\n')
				l.write(hit+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t'+str(E)+'\t'+str(S)+'\t'+QUERY+'\n')
				print hit+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t'+str(E)+'\t'+str(S)+'\t'+QUERY
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
	
