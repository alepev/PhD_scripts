#!/usr/bin/python

import sys

if len(sys.argv)!=5:
	print 'usage: <blast output> <fasta file> <score threshold> <e-value threshold>'
else:
	blast= sys.argv[1]
	fasta= sys.argv[2]
	S=int(sys.argv[3])
	E=float(sys.argv[4])
	fastaout= '.'.join(blast.split('/')[-1].split('.')[:-1])+'.out.fa'
	logfile= '.'.join(blast.split('/')[-1].split('.')[:-1])+'.out.log'

	#load blast results

	HITS= {}
	READ= 0
	for line in open(blast).readlines():
		if READ>0:
			if line.split()==[]: #initial empty line ignored, final causes READ==0
				READ-=1
			else:
				x= line.split('|')
				if len(x)==1:
					ID= line.split()[0]
				else:
					ID= '|'.join(x[:-1])
				block= line[70:].lstrip().split()
				score= int(block[0])
				if block[1][0]=='e':
					evalue= float('1'+block[1])
				else:
					evalue= float(block[1])
				if score >= S and evalue <= E:
					if HITS.get(ID,0)==0:
						HITS[ID]=(score,evalue)
					elif HITS[ID][1]>evalue or (HITS[ID][1]==evalue and HITS[ID][0]<score):
						HITS[ID]=(score,evalue)
		elif len(line.split('Sequences producing significant alignments:'))>1:
			READ=2

	#parse fasta seqs

	f=open(fastaout,'w')
	APPEND=True
	for line in open(fasta).readlines():
		x= line.split()
		if line[0]=='>':
			locus=x[0][1:]
			if len(locus)>=65:
				locus=locus[:65]
			if len(locus.split('|'))>1:
				locus='|'.join(locus.split('|')[:-1])
			if HITS.get(locus,0)!=0:
				APPEND=True
				f.write(line)
			else:
				APPEND=False				
		elif APPEND:
			f.write(line)
	f.close()

	#output summary

	hits= HITS.keys()
	hits.sort()
	f= open(logfile,'w')
	f.write('#hit\tscore\te-value\n')
	for h in hits:
		values = '\t%d\t%.2e' % (HITS[h][0], HITS[h][1])
		f.write(h+values+'\n')
	f.close() 

	
