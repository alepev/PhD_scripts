#!/usr/bin/python

#v3: FASTA only for best hits!
#IMPORTANT NOTE: ADDED SUPPORT FOR SPLICED SEQUENCES (curly brackets), BUT NOT SURE IT WORKS!

import sys
import subprocess as sub
sub.PIPE=1

AA= {}
AA['Ala']='A'
AA['Arg']='R'
AA['Asn']='N'
AA['Asp']='D'
AA['Cys']='C'
AA['Glu']='E'
AA['Gln']='Q'
AA['Gly']='G'
AA['His']='H'
AA['Ile']='I'
AA['Leu']='L'
AA['Lys']='K'
AA['Met']='M'
AA['Phe']='F'
AA['Pro']='P'
AA['Ser']='S'
AA['Thr']='T'
AA['Trp']='W'
AA['Tyr']='Y'
AA['Val']='V'
AA['Unk']='X'

if len(sys.argv) not in [2,3]:
	print
	print 'usage: <EXONERATE output> <score threshold> (OPTIONAL: default ALL hits)'
	print
	print '	OUTPUT:'
	print '	 - hits sequence limits file (BEST + other hits commented)'
	print '	 - hits aligned regions fasta file, DNA (BEST only)'
	print '	 - hits aligned regions fasta file, PROTEIN (BEST only)'	
	print 

#ABOUT EXONERATE INDEXING:
#limits reported next to sequences are inclusive indexes of first-last residues in a system starting with index 1;
#the output preserves the original indexing system,
#so sequence extraction operations require a conversion from start:end (inclusive/inclusive) -> start-1:end (inclusive/exclusive)
#ex. sequence 1:4 (indicating residues 1-2-3-4) -> string[0:4] (indicating residues 0-1-2-3, where 0 corresponds to previos 1)

else:
	infile= sys.argv[1]
	if len(sys.argv)==3:
		S= int(sys.argv[2])
		ALL=False
	else:
		ALL=True
		
	ntout= infile+'.nt.fa'
	aaout= infile+'.aa.fa'
	limitout= infile+'.limits.txt'

	#extract Exonerate results
	HITS={} 
	RESULTS= '\n'.join(open(infile).readlines()).split('Query:')
	if len(RESULTS)>1:
		for q in RESULTS[1:]:
			QUERY= q.split()[0]
			#target sequence features
			HIT= q.split('Target: ')[1].split('\n')[0].split(':[revcomp]')
			if len(HIT)>1:
				STRAND='-'
			else:
				STRAND='+'
			HIT=HIT[0].split('len=')[0]
			SCORE= int(q.split('Raw score: ')[1].split()[0])
			if not ALL and SCORE<S:
				continue
			#THE FOLLOWING IS DISMISSED, HAVE TO CHECK HOW INDEXING IS HANDLED FOR [-]!
			#x= q.split('Target range: ')[1].split(' -> ') 
			#START= int(x[0].split()[-1])
			#STOP= int(x[1].split()[0])
			x= q.split('Target range: ')[1].split('vulgar:')[0].split('\n')[1:]
			aaseq=[]
			ntseq=[]
			for line in x:
				if len(line.split())!=0:
					x1= line.lstrip()[0]
					if x1.isalpha() or x1 in ['*','-','{','}','>','<'] or (x1=='#' and len(line.split('|'))==1): 
						aaseq.append(line.lstrip().rstrip())
					elif x1.isdigit():
						ntseq+=[line.lstrip()]
			#nt seq
			NTseq=''
			START=0
			for i in range(len(ntseq)):
				if (i+1)%2 == 0: #only even lines are the ones to keep, the others contain query AA seq
					NTseq+=ntseq[i].split()[2]
					if START==0:
						START=int(ntseq[i].split()[0]) #only first line
					STOP=int(ntseq[i].split()[4]) #always, so value in last iteration is the one kept					
			NTseq=NTseq.replace('-','') #remove gaps
			#aa seq
			aaseq=''.join(aaseq).replace('-','') #remove gaps
			i=0
			aalen=0
			AAseq=''
			while i <(len(aaseq)-2):
				if aaseq[i].isupper():
					AAseq+=AA[aaseq[i:i+3]]
					i+=3
					aalen+=1
				elif aaseq[i]=='*':
					AAseq+='*'
					i+=3
				else:
					i+=1 #in case of '#'
			if START>STOP:
				START,STOP=STOP,START
			HITS[HIT]= HITS.get(HIT,[])+[(SCORE,START,STOP,STRAND,QUERY,NTseq,AAseq,aalen)]

	#output summary
	if len(HITS.keys())>0:
		n= open(ntout,'w')
		a= open(aaout,'w')
		l= open(limitout,'w')
		l.write('#ID\tSTART\tSTOP\tSTRAND\tNTlen\tAAlen\tEx.SCORE\tQUERY\n')
		print '#ID\tSTART\tSTOP\tSTRAND\tNTlen\tAAlen\tEx.SCORE\tQUERY'
		SORTED=[]
		for HIT in HITS.keys():
			scores=HITS[HIT]
			scores.sort(reverse=True)
			HITS[HIT]=scores
			SORTED+= [(scores[0][0],HIT)] #best S for hit
		SORTED.sort(reverse=True)
		for sorted in SORTED:
			HIT=sorted[1]
			comment=False #all scores except 1st (BEST) for a given hit are commented in limits file
			for score in HITS[HIT]:
				SCORE=score[0]
				START= score[1]
				STOP= score[2]
				STRAND= score[3]
				QUERY= score[4]
				NTseq= score[5]
				AAseq= score[6]
				aalen= score[7]
				if comment:
					l.write('#'+HIT+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t('+str(aalen)+')\t'+str(SCORE)+'\t'+QUERY+'\n')
					print '#'+HIT+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t('+str(aalen)+')\t'+str(SCORE)+'\t'+QUERY
				else:
					l.write(HIT+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t('+str(aalen)+')\t'+str(SCORE)+'\t'+QUERY+' (BEST)\n')
					print HIT+'\t'+str(START)+'\t'+str(STOP)+'\t'+STRAND+'\t('+str(STOP-START+1)+')\t('+str(aalen)+')\t'+str(SCORE)+'\t'+QUERY+' (BEST)'
					n.write('>'+HIT+' | '+str(START)+':'+str(STOP)+' ['+STRAND+'] ('+str(STOP-START+1)+'nt)- Exonerate vs '+QUERY+' (S'+str(SCORE)+')\n')
					n.write(NTseq+'\n')
					a.write('>'+HIT+' | '+str(START)+':'+str(STOP)+' ['+STRAND+'] ('+str(aalen)+'aa)- Exonerate vs '+QUERY+' (S'+str(SCORE)+')\n')
					a.write(AAseq+'\n')
				comment=True
				#print NTseq
				#print AAseq
		n.close()
		a.close()
		l.close()

		print
		print '*********************************************************'
		print 'check output files for results:'
		print ' - aligned sequences limits:',limitout
		print ' - aligned sequences fasta, DNA:', ntout
		print ' - aligned sequences fasta, PROTEIN:', aaout
		print
	else:
		print
		print '*** no Exonerate hit retrieved! ***'
		print
	
