#!/usr/bin/python

import sys
import copy

#v2.0: uses partial names from tabular output, NOT alignment as in previous version because they are not shown for all hits!

EVAL=0.000000001

if len(sys.argv) not in [3,4]:
	print
	print 'usage: <BLASTP+ output> <PROTEOME fasta> <SCORE threshold> (OPTIONAL: default ALL hits)'
	print 'NOTE: assumes blast+ output from search with -DB OPTION, standard format'
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
		QUERY= ''.join(query.split('Length=')[0].split('\n'))
		y= query.split('Identities =') #-db output option
		if len(y)==1:
			print '*** no hits found for',QUERY,'***'
		else:
			print y[0]
			hits=y[0].split('Value')[1].split('\n>')[0].split('\n')
			for hit in hits:
				if hit!='':
					HIT= hit[2:65].rstrip()	
					print HIT[68:]			
					S= float(hit[68:].split()[0])
					E= float(hit[68:].split()[1])
					if not ALL and S<score or E>EVAL:
						#print hit[68:]
						#print HIT,S,E
						continue #HIT SCORE BELOW THRESHOLD
					if HITS.get(HIT,0)==0:
						HITS[HIT]=(S,E,QUERY)
					elif S > HITS[HIT][0]:
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
			if len(x)>0 and x[0]=='>':
				locus=x[1:64]
				if HITS.get(locus,0)!=0:
					APPEND=True
					Nfound[locus]=1
					f.write(x+'\n')
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

