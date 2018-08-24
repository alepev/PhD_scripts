#!/usr/bin/python

#WHAT IS THIS PROGRAM FOR:
#given a list of loci, extacts sequence context of the main ORF START (ATG)
#(how many nucleotides before and after can be selected below in the control panel)
#uses coordinates from gff file combined with transcript (cDNA) sequence

#================================================

#CONTROL PANEL

ATGup=5		#INCLUSIVE (ex. 5 -> -5 to -1)
ATGdw=8 	#INCLUSIVE (ex. 5 -> +1 to +5; NOTE: includes ATG, i.e. +5 = ATG + 2 nts, not ATG + 5 nts!)

modelfile=	'TAIR10_representative_gene_models.txt'
gff3file=	'TAIR10_GFF3_genes.gff'
seqfile=	'TAIR10_seq_20101214_updated.fa' #includes introns
#NOTE: tried with cDNA first, which caused coordinate shift for many transcripts, i.e. invalid start codon

#================================================

#load target list

import sys

if len(sys.argv)!=2:
	print
	print 'usage: <gene list (AGI loci)>'
	print
else:
	infile= sys.argv[1]
	out1='.'.join(infile.split('.')[:-1])+'.-'+str(ATGup)+'+'+str(ATGdw)+'.seqs.out'	#context sequence file
	out2='.'.join(infile.split('.')[:-1])+'.-'+str(ATGup)+'+'+str(ATGdw)+'.freq.out'	#frequencies file
	
#================================================

#load input gene models
	
	MOD={}
	for line in open(modelfile,'rU').readlines():
			x=line.split()
			if len(x)!=0 and x[0][:2]=='AT':
				locus=x[0].split('.')
				MOD[locus[0]]=locus[1]
	
	NOTFOUND={}
	LOCI={}
	for line in open(infile,'rU').readlines():
		x=line.split()[0]
		if len(x)>0 and x[0]!='#':
			if x[0]=='>':
				x=x[1:]
			if x[:2]=='AT':
				locus=x.split('|')[0].split('/')[0] #remove descriptors, add here other separators if needed
				if len(locus.split('.'))>1: #gene model provided
					LOCI[locus]=1
				elif MOD.get(locus,0)!=0: #representative gene model automatically assigned
					LOCI[locus+'.'+MOD[locus]]=1
				else:
					NOTFOUND[locus]=1

	Nloci=len(LOCI.keys())
	Nnomodel= len(NOTFOUND.keys())

#================================================

#extract coordinates

	START={}		#ATG position on transcript
	SHORT={}		#not possible to extract consensus because of too short upstream/downstream sequence around ATG	
	
	READ=False
	for line in open(gff3file,'rU').readlines()[1:]:
		x=line.split()
		feat=x[2]
		if feat=='mRNA':	
			locus=x[8].split(';')[0][3:]
			if LOCI.get(locus,0)!=0:
				strand=x[6]
				if strand=='+':
					start1=int(x[3])
					stop1=int(x[4])
				else:
					stop1=int(x[3])
					start1=int(x[4])
				L=abs(stop1-start1)+1
				READ=True
		elif feat=='CDS' and READ:
			#locus=x[8].split(',')[0][7:]
			if strand=='+':
				start2=int(x[3])
				stop2=int(x[4])
			else:
				start2=int(x[4])
				stop2=int(x[3])
			relativeStart = abs(start2-start1)
			if relativeStart >= ATGup and (L - relativeStart +1) >= ATGdw:
				START[locus]=relativeStart
			else:
				SHORT[locus]=1
			READ=False

	Ncoding=len(START.keys())
	Nshort=len(SHORT.keys())

#================================================

#extract ATG context from transcript sequence

	SEQ={}
	BADSTART={}
			
	locus=''
	seq=''
	READ=False
	for line in open(seqfile,'rU').readlines():
			if line[0]=='>':
				if locus!='' and READ: #if new ID, check previous sequence
					if seq[START[locus]:START[locus]+3]=='ATG':
						SEQ[locus]=seq[START[locus]-ATGup:START[locus]+ATGdw]
						#print seq[START[locus]-ATGup:START[locus]+ATGdw]
					else:
						BADSTART[locus]=1
				locus=line.split()[0][1:] #save new sequence ID
				seq=''
				if START.get(locus,0)!=0:
					READ=True
				else:
					READ=False			
			elif READ: #if sequence and corresponding ID was OK, read it
				if seq=='':
					seq=line.split()[0].upper()
				else:
					seq+=line.split()[0].upper() #add to existing
						
	#last sequence in the file
	if locus!='' and READ:
		if seq[START[locus]:START[locus]+3]=='ATG':
			SEQ[locus]=seq[START[locus]-ATGup:START[locus]+ATGdw]
			#print seq[START[locus]-ATGup:START[locus]+ATGdw]
		else:
			BADSTART[locus]=1		
	
	Nseqs=len(SEQ.keys())
	Nbad=len(BADSTART.keys())

#================================================

#compute frequencies at each position of ATG context

	loci=SEQ.keys()
	loci.sort()

	COUNTS={}
	for i in range(ATGup+ATGdw):
		COUNTS[i]={'A':0,'C':0,'G':0,'T':0,'N':0} #(filled below)

	f=open(out1,'w')
	#f.write('#locus\tcontext_sequence')
	for locus in loci:
		#f.write('\n'+locus+'\t'+SEQ[locus])
		f.write('>'+locus+'\n'+SEQ[locus]+'\n')
		for i in range(len(SEQ[locus])):
			nt=SEQ[locus][i]
			COUNTS[i][nt]=COUNTS[i][nt]+1
	#		if COUNTS[i].has_key(nt):
	#			COUNTS[i][nt]=COUNTS[i][nt]+1
	#		else:
	#			COUNTS[i]['N']=COUNTS[i]['N']+1 #in theory N already included, but in case of other characters
	f.close()

	f=open(out2,'w')
	f.write('\n#FILTERING')
	f.write('\nTOTAL INPUT LOCI:\t'+str(Nloci+Nnomodel))
	f.write('\n-without TAIR model:\t'+str(Nnomodel))
	f.write('\n=with TAIR model:\t'+str(Nloci))
	f.write('\n-not protein-coding:\t'+str(Nloci-(Nseqs+Nshort)))
	f.write('\n=protein-coding:\t'+str(Nseqs+Nshort))
	f.write('\n-too short to extract context:\t'+str(Nshort))
	f.write('\n=sufficient length:\t'+str(Ncoding))
	f.write('\n-invalid start codon:\t'+str(Nbad))
	f.write('\n=VALID ATG SEQUENCE CONTEXT:\t'+str(Nseqs))
	f.write('\n\n#NT COUNTS\n')
	for i in range(-ATGup,ATGdw,1):
		if i<0:
			f.write('\t'+str(i))
		else:
			f.write('\t+'+str(i+1))
	for nt in ['A','C','G','T','N']:
		f.write('\n'+nt)
		for i in range(ATGup+ATGdw):
			f.write('\t'+str(COUNTS[i][nt]))
	f.write('\n\n#NT FREQ.s\n')
	for i in range(-ATGup,ATGdw,1):
		if i<0:
			f.write('\t'+str(i))
		else:
			f.write('\t+'+str(i+1))
	for nt in ['A','C','G','T','N']:
		f.write('\n'+nt)
		for i in range(ATGup+ATGdw):
			f.write('\t'+str(COUNTS[i][nt]*100.0/Nseqs))
	f.close()
	

