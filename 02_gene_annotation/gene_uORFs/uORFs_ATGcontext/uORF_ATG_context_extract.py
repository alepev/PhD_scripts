#!/usr/bin/python

#WHAT IS THIS PROGRAM FOR:
#given a list of loci, extracts context of uORF START (ATG) for:
#- ALL uORFs
#- only strongest/longest/most overlapping for each gene (consider also frame?)
#uses data from supplementary table in von Arnim et al. 2014

#================================================

#CONTROL PANEL

modelfile=	'TAIR10_representative_gene_models.txt'
gff3file=	'TAIR10_GFF3_genes.gff'
uORFfile= 	'von_Arnim_2014_TAIR10_uORFs.txt'

#================================================

#load target list

import sys

if len(sys.argv)!=2:
	print
	print 'usage: <gene list (AGI loci)>'
	print
else:
	infile= sys.argv[1]
	out1a='.'.join(infile.split('.')[:-1])+'.ATG_ORFall.fa'	#context sequence file, ALL uORFs
	out1b='.'.join(infile.split('.')[:-1])+'.ATG_1ORF.fa'		#context sequence file, best uORF per gene
	out2='.'.join(infile.split('.')[:-1])+'.ATG_ORF.log'		#frequencies file
	
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

#extract ATG context from data table

	SEQ={}	#SEQ of all uORFs
	BEST={}	#best uORF per gene
	BADorf=0
	
	for line in open(uORFfile,'rU').readlines()[1:]:
		x=line.split()
		locus=x[0]
		uORF=x[1]
		if LOCI.get(locus,0)==0:
			continue
		elif uORF=='-':
			LOCI[locus]=-1 #mapped, no uORFs
		else:
			ATGseq=x[7]
			if len(ATGseq)<13:
				BADorf+=1
				if LOCI[locus]!=2:
					LOCI[locus]=-2	#no valid legth ATG context			
			else:
				LOCI[locus]=2  #mapped, with at least 1 valid uORF(s)
				frame=int(x[4])
				ovlength=int(x[6])
				strength=int(x[8])
				length=int(x[9])
				SEQ[uORF]=ATGseq
				if BEST.get(locus,0)==0:
					BEST[locus]=(uORF,strength,length,ovlength,frame)
				else:
					s,l,ovl,f=BEST[locus][1],BEST[locus][2],BEST[locus][3],BEST[locus][4]
					if strength>s:
						BEST[locus]=(uORF,strength,length,ovlength,frame)
					elif strength == s:
						if length > l:
							BEST[locus]=(uORF,strength,length,ovlength,frame)
						elif length == l:
							if ovlength > ovl:
								BEST[locus]=(uORF,strength,length,ovlength,frame)
							elif ovlength == ovl:
								if frame < f:
									BEST[locus]=(uORF,strength,length,ovlength,frame)

									
#================================================

#compute frequencies at each position of ATG context
										
	COUNTS={}
	COUNTSBEST={}
	for i in range(13):
		COUNTS[i]={'A':0,'C':0,'G':0,'T':0,'N':0} 
		COUNTSBEST[i]={'A':0,'C':0,'G':0,'T':0,'N':0}

	uorfs=SEQ.keys()
	uorfs.sort()
	fa=open(out1a,'w')
	fb=open(out1b,'w')
	for u in uorfs:
			fa.write('>'+u+'\n'+SEQ[u]+'\n')
			for i in range(13):
				nt=SEQ[u][i]
				COUNTS[i][nt]=COUNTS[i][nt]+1
			if BEST['.'.join(u.split('.')[:-1])][0]==u:
				fb.write('>'+u+'\n'+SEQ[u]+'\n')
				for i in range(13):
					nt=SEQ[u][i]
					COUNTSBEST[i][nt]=COUNTSBEST[i][nt]+1
	fa.close()
	fb.close()
	
	NorfTOTAL=len(uorfs)
	Norfloci=len(BEST.keys())
	Nnomap= len(filter(lambda x:LOCI[x]==1,LOCI))
	Nnoorf= len(filter(lambda x:LOCI[x]==-1,LOCI))
	Nshort= len(filter(lambda x:LOCI[x]==-2,LOCI))
	
#================================================

#compute frequencies at each position of ATG context

	f=open(out2,'w')
	f.write('\n#FILTERING')
	f.write('\nTOTAL INPUT LOCI:\t'+str(Nloci+Nnomodel))
	f.write('\n-without TAIR model:\t'+str(Nnomodel))
	f.write('\n=with TAIR model:\t'+str(Nloci))
	f.write('\n-unmapped:\t'+str(Nnomap))
	f.write('\n-no uORFs:\t'+str(Nnoorf))
	f.write('\n-no valid length uORFs:\t'+str(Nshort))
	f.write('\n=WITH uORF(s):\t'+str(Norfloci))
	f.write('\n\nTOTAL uORFs:\t'+str(NorfTOTAL))
	f.write('\n\n#NT COUNTS (ALL)\n')
	for i in range(-5,8,1):
		if i<0:
			f.write('\t'+str(i))
		else:
			f.write('\t+'+str(i+1))
	for nt in ['A','C','G','T','N']:
		f.write('\n'+nt)
		for i in range(13):
			f.write('\t'+str(COUNTS[i][nt]))
	f.write('\n\n#NT FREQ.s (ALL)\n')
	for i in range(-5,8,1):
		if i<0:
			f.write('\t'+str(i))
		else:
			f.write('\t+'+str(i+1))
	for nt in ['A','C','G','T','N']:
		f.write('\n'+nt)
		for i in range(13):
			f.write('\t'+str(COUNTS[i][nt]*100.0/NorfTOTAL))
	f.write('\n\n#NT COUNTS (BEST)\n')
	for i in range(-5,8,1):
		if i<0:
			f.write('\t'+str(i))
		else:
			f.write('\t+'+str(i+1))
	for nt in ['A','C','G','T','N']:
		f.write('\n'+nt)
		for i in range(13):
			f.write('\t'+str(COUNTSBEST[i][nt]))
	f.write('\n\n#NT FREQ.s (BEST)\n')
	for i in range(-5,8,1):
		if i<0:
			f.write('\t'+str(i))
		else:
			f.write('\t+'+str(i+1))
	for nt in ['A','C','G','T','N']:
		f.write('\n'+nt)
		for i in range(13):
			f.write('\t'+str(COUNTSBEST[i][nt]*100.0/Norfloci))
	f.write('\n\n(BEST ORF: strength > length > overlap_length > frame)')
	f.close()

	

