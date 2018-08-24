#!/usr/bin/python

import sys
import copy

EVAL=0.001

if len(sys.argv) not in [3,4]:
	print
	print 'usage: <BLASTP+ output> <PROTEOME fasta> <SCORE threshold> (OPTIONAL: default ALL hits)'
	print 'NOTE:'
	print '- assumes blast+ output from search with -DB OPTION, standard format'
	print '- uses full names from blast alignments to match fasta file'
	print '(slower, but ensures original name matching, which is frequently pruned in tabular summary)'
	print
	print '	OUTPUT:'
	print '	 - FASTA containing hits sequences, selected according to SCORE threshold (if provided)'
	print	'		and E-value threshold of',EVAL,' (see script to modify)'
	print 

else:
	blast= sys.argv[1]
	fasta= sys.argv[2]
	if len(sys.argv)==4:
		score= float(sys.argv[3])
		ALL=False
		fastaout= '.'.join(blast.split('.')[:-1])+'.S'+str(int(score))+'.out.fa'
		summary= '.'.join(blast.split('.')[:-1])+'.S'+str(int(score))+'.out.log'
	else:
		ALL=True
		fastaout= '.'.join(blast.split('.')[:-1])+'.out.fa'
		summary= '.'.join(blast.split('.')[:-1])+'.out.log'
	
	#extract blast results
	HITS= {}
	x= '\n'.join(open(blast).readlines()).split('Query= ')
	for query in x[1:]:
		QUERY= query.split('\n')[0].rstrip()
		y= query.split('> ') #-db output option
		if len(y)==1:
			print '*** no hits found for',QUERY,'***'
			continue #NO HIT FOR THAT QUERY
		for subj in y[1:]: #many hits can be reported for same subject when long target sequences 
			HIT= ''.join(subj.split('Length=')[0].split('\n')) #merges long names on multiple lines
			y2= subj.split('Score = ')
			for hit in y2[1:]:
				S= float(hit.split(' bits')[0])
				E= float(hit.split('Expect')[1].split(',')[0].split()[-1]) #normally "Expect =", however when multiple hits in -db option -> "Expect(N) = "
				if not ALL and S<score or E > EVAL:
					continue #HIT SCORE BELOW THRESHOLD
				if HITS.get(HIT,0)==0:
					HITS[HIT]=(S,E,QUERY)
				else:
					s1=HITS[HIT][0]
					#e1=HITS[HIT][1]
					if S>s1:
						HITS[HIT]=(S,E,QUERY) 
						
	if HITS!={}:		
		OUTPUT=[]
		for HIT in HITS.keys():
			OUTPUT+=[(-HITS[HIT][0],HITS[HIT][1],HIT,HITS[HIT][2])] #-score to allow sort from high to low, as for e-value
		OUTPUT.sort()

		#output summary
		s= open(summary,'w')
		s.write('#QUERY\tHIT\tSCORE\tE-VALUE\n')
		print
		print '#QUERY\tHIT\tSCORE\tE-VALUE'
		for OUT in OUTPUT:
			S,E,HIT,QUERY=-OUT[0],OUT[1],OUT[2],OUT[3]
			s.write(QUERY+'\t'+HIT+'\t'+str(S)+'\t'+str(E)+'\n')
			print QUERY+'\t'+HIT+'\t'+str(S)+'\t'+str(E)
		s.close()
		
		#parse fasta seqs
		
		Nfound={}
		f=open(fastaout,'w')
		APPEND=True
		for line in open(fasta).readlines():
			x= line.split('\n')[0]
			if x[0]=='>':
				locus=x[1:]
				if HITS.get(locus,0)!=0:
					APPEND=True
					Nfound[locus]=1
					f.write('>'+locus+'\n')
				else:
					APPEND=False				
			elif APPEND:
				f.write(line)
		f.close()
		
		if len(Nfound)!=len(HITS):
			print '\n***WARNING: the following',len(HITS)-len(Nfound),'hits could not be parsed!***'
			print '(missing in provided fasta file / mispelled locus name)'
			for locus in HITS.keys():
				if Nfound.get(locus,0)==0:
					print locus
	
	else:
		print '\n***NO HITS FOUND***\n'

