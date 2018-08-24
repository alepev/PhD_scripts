#!/usr/bin/python

import sys

if len(sys.argv) not in [3,5]:
	print '\nusage: <blast output> <fasta file> <score threshold> <e-value threshold>'
	print '(FASTA PARSING ONLY: omit score and e-value threshold (no log file generated))'
	print '\nNOTE: assumes blast+ output with lines containing only subject, score, e-value'
	print '(as generated by -outfmt "6 sseqid score evalue" in blast+ settings)\n'
else:
	blast= sys.argv[1]
	fasta= sys.argv[2]
	if len(sys.argv)==5:
		S=int(sys.argv[3])
		E=float(sys.argv[4])
		fastaout= '.'.join(blast.split('/')[-1].split('.')[:-1])+'.out'
		logfile= '.'.join(blast.split('/')[-1].split('.')[:-1])+'.out'
		if S!=0:
			fastaout= fastaout+str(S)
			logfile= logfile+str(S)
		if E!=1:
			xx='%.2e' % E
			fastaout= fastaout+xx[-4:]
			logfile= logfile+xx[-4:]
		fastaout= fastaout+'.fa'
		logfile= logfile+'.log'
	else:
		S=0
		E=1
		fastaout= '.'.join(blast.split('/')[-1].split('.')[:-1])+'.out.fa'

	#load blast results

	HITS= {}
	for line in open(blast).readlines():
		x= line.split('\t')
		ID=x[0]
		#if len(ID.split('|'))==1:
		#	ID= ID.split()[0]
		#else:
		#	ID= '|'.join(ID.split('|')[:-1])
		if S!=0 and E!=1:
			score=int(x[1])
			if x[2][0]=='e':
				evalue= float('1'+x[2][:-1])	#removes final '\n'
			else:
				evalue= float(x[2][:-1])	#removes final '\n'
			if score >= S and evalue <= E:
				if HITS.get(ID,0)==0:
					HITS[ID]=(score,evalue)
				elif HITS[ID][1]>evalue or (HITS[ID][1]==evalue and HITS[ID][0]<score):
					HITS[ID]=(score,evalue)
		else:
			HITS[ID]=1

	#parse fasta seqs

	f=open(fastaout,'w')
	APPEND=True
	for line in open(fasta).readlines():
		x= line.split()
		if line[0]=='>':
			locus=x[0][1:]
			#if len(locus.split('|'))>1:
			#	locus='|'.join(locus.split('|')[:-1])
			if HITS.get(locus,0)!=0:
				APPEND=True
				f.write('>'+locus+'\n')
			else:
				APPEND=False				
		elif APPEND:
			f.write(line)
	f.close()

	#output summary

	if S!=0 and E!=1:
		hits= HITS.keys()
		hits.sort()
		f= open(logfile,'w')
		f.write('#hit\tscore\te-value\n')
		for h in hits:
			values = '\t%d\t%.2e' % (HITS[h][0], HITS[h][1])
			f.write(h+values+'\n')
		f.close() 

	
