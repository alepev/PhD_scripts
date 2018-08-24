#!/usr/bin/python

import sys

if len(sys.argv) != 2:
	print 'usage: <FASTA file> from NCBI'
	print 'moves species name at the beginning (if not found, indicated by  ??)'
else:
	fasta= sys.argv[1]
	out= fasta+'.map'
	log= fasta+'.log'

	MAP={}
	BAD=[]
	unmapped=0
	x= open(fasta).readlines()
	f=open(out,'w')
	for i in range(len(x)-1):
		if x[i][0]=='>':
			xs= x[i].split()
			if xs[-1].isalpha() and xs[-1].isupper() and len(xs[-1])>2:
				header= ' '.join(xs[:-1])
				seq= 	'\n'+xs[-1]
			else:
				header= x[i][:-1] #excludes final '\n'
				seq='\n'
				if not (x[i+1]!='\n' and x[i+1][0]!='>'):
					BAD.append(header)
					continue
			if header.split()[-1][-1]==']':
				code= '_'.join(header.split('[')[-1].split(']')[0].split())
				MAP[code]=1
				f.write('>'+code+'_'+header[1:header.rfind('[')]) 
			elif len(header.split('|')[-1].split()[0].split('_'))==2:
				code= header.split('|')[-1].split()[0].split('_')[1]
				f.write('>'+code+'_'+header[1:]) 
				MAP[code]=1
			else:
				f.write('>??_'+header[1:]) 
				unmapped+=1
			f.write(seq)
		elif x[i]!='\n': #removes empty lines
			f.write(x[i])
	f.close()

	m= MAP.keys()
	m.sort()
	f=open(log,'w')
	if unmapped>0:
		f.write('#unmapped sequences (see >?? in fasta output):\t'+str(unmapped)+'\n')
	if len(BAD)>0:
		f.write('#seqIDs with no corresponding sequence:\t'+str(len(BAD))+'\n')
		for b in BAD:
			f.write(b+'\n')
	f.write('#MAPPED SEQUENCES\n')
	for code in m:
		f.write(code+'\n')
	f.close()
