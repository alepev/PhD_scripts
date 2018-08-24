#!/usr/bin/python

import sys

CUT=50 #generate counts only for the first X codons of the CDS; IGNORED IF == 0

CDSfile= 'TAIR10_cds_20101214_updated.fa'
modelfile='TAIR10_representative_gene_models.txt'

if len(sys.argv)>2:
	print
	print 'usage: COUNT CODONS <ID list (OPTIONAL)>'
	print
	print 'NOTE:'
	print '- start (1st Met == ATG) and stop codons (TAA,TAG,TGA) are ignored in the statistics'
	print '- bad codons (containing a character other than A,C,G,T) are also ignored'
	print '- CDSs with invalid length (!= 3N) are trimmed at the 3-end'
	print '- CDSs with invalid start (!= ATG) are completely excluded from the analysis'

	if CUT>0:
		print '- only first',CUT,'CDS codons will be used for statistics (see code to change this value)' 
	print
else:

	#load representative gene models

	MOD={}
	for line in open(modelfile,'rU').readlines():
			x=line.split()
			if len(x)!=0 and x[0][:2]=='AT':
				locus=x[0].split('.')
				MOD[locus[0]]=locus[1]
				
	#load loci list

	LOCI={}
	if len(sys.argv)==2:
		listfile=sys.argv[1]
		NOTFOUND=[]
		for line in open(listfile,'rU').readlines():
			x=line.split()[0]
			if len(x)>0 and x[0]!='#' and x[:2]=='AT':
				locus=x.split('|')[0].split('/')[0] #remove descriptors, add here other separators if needed
				if len(locus.split('.'))>1: #gene model provided
					LOCI[locus]=1
				elif MOD.get(locus,0)!=0: #representative gene model automatically assigned
					LOCI[locus+'.'+MOD[locus]]=1
				else:
					NOTFOUND.append(locus)
	
		if NOTFOUND!=[]:
			print
			print '***WARNING: the following IDs could not be mapped to a gene model (excluded from parsing):'
			for locus in NOTFOUND:
				print locus
			print
	
	#####################################################
	#STANDARD TRANSLATION CODE
	#####################################################
	CODE= {}
	CODE['ATG']='Met' #translation start
	CODE['TGG']='Trp'
	for triplet in ['TTT','TTC']:
		CODE[triplet]='Phe'
	for triplet in ['TAT','TAC']:
		CODE[triplet]='Tyr'
	for triplet in ['TGT','TGC']:
		CODE[triplet]='Cys'
	for triplet in ['AAT','AAC']:
		CODE[triplet]='Asn'
	for triplet in ['CAA','CAG']:
		CODE[triplet]='Gln'
	for triplet in ['GAT','GAC']:
		CODE[triplet]='Asp'
	for triplet in ['GAA','GAG']:
		CODE[triplet]='Glu'
	for triplet in ['AAA','AAG']:
		CODE[triplet]='Lys'
	for triplet in ['CAT','CAC']:
		CODE[triplet]='His'
	for triplet in ['ATT','ATC','ATA']:
		CODE[triplet]='Ile'
	for triplet in ['GTT','GTC','GTA','GTG']:
		CODE[triplet]='Val'
	for triplet in ['ACT','ACC','ACA','ACG']:
		CODE[triplet]='Thr'
	for triplet in ['GCT','GCC','GCA','GCG']:
		CODE[triplet]='Ala'
	for triplet in ['CCT','CCC','CCA','CCG']:
		CODE[triplet]='Pro'
	for triplet in ['GGT','GGC','GGA','GGG']:
		CODE[triplet]='Gly'
	for triplet in ['TTA','TTG','CTT','CTC','CTA','CTG']:
		CODE[triplet]='Leu'
	for triplet in ['TCT','TCC','TCA','TCG','AGT','AGC']:
		CODE[triplet]='Ser'
	for triplet in ['CGT','CGC','CGA','CGG','AGA','AGG']:
		CODE[triplet]='Arg'
	for triplet in ['TAA','TAG','TGA']:
		CODE[triplet]='STOP' #translation stop
	
	#scan CDS file for aa/triplet counts

	TRIPLET={}
	AA={}
	NAME=''
	seq=''
	BADstart={}
	BADlen={}
	COUNT=False
	TOTcod=0
	BADcod=0
	Nseqs=0
	GC1=0 #1st codon position
	GC2=0	#2nd codon position
	GC3=0	#3rd codon position
	for line in open(CDSfile,'rU').readlines():
			if line[0]=='>':
				if NAME!='' and COUNT:
					if len(seq)%3!=0:
						BADlen[NAME]=1
						seq=seq[:len(seq)-len(seq)%3]
						if CUT==0:
							MAX=len(seq)
						else:
							MAX=min(len(seq),CUT*3)
						for i in range(0,MAX,3):
							triplet=seq[i:i+3]
							if CODE.get(triplet,0)!=0:
								if CODE.get(triplet)!='STOP': #stop codons are not relevant (1 per seq)
									TRIPLET[triplet]=TRIPLET.get(triplet,0)+1
									AA[CODE[triplet]]=AA.get(CODE[triplet],0)+1
									TOTcod+=1
									if triplet[0] in ['G','C'] :
										GC1+=1
									if triplet[1] in ['G','C']:
										GC2+=1
									if triplet[0] in ['G','C']:
										GC3+=1
							else:
								BADcod+=1
				NAME=line.split()[0][1:]
				lc=0
				seq=''
				if LOCI=={} or LOCI.get(NAME,0)!=0:
					Nseqs+=1
					COUNT=True
				else:
					COUNT=False
			elif COUNT:
				if seq=='':
					start=line[:3].upper()
					if start!='ATG':
						 BADstart[NAME]=1
						 COUNT=False
					else:
						seq=line.split()[0][3:].upper() #1st ATG is not relevant (1 per seq)
				else:
						seq+=line.split()[0].upper()
	if NAME!='' and COUNT:
		if len(seq)%3!=0:
			BADlen[NAME]=1
	else:
		if CUT==0:
			MAX=len(seq)
		else:
			MAX=min(len(seq),CUT*3)
		for i in range(0,MAX,3):
			triplet=seq[i:i+3]
			if CODE.get(triplet,0)!=0:
				if CODE.get(triplet)!='STOP': #stop codons are not relevant (1 per seq)
					TRIPLET[triplet]=TRIPLET.get(triplet,0)+1
					AA[CODE[triplet]]=AA.get(CODE[triplet],0)+1
					TOTcod+=1
					if triplet[0] in ['G','C'] :
						GC1+=1
					if triplet[1] in ['G','C']:
						GC2+=1
					if triplet[0] in ['G','C']:
						GC3+=1
			else:
				BADcod+=1
				
	#output statistics
	
	AAuniq={}
	for aa in CODE.values():
		AAuniq[aa]=1
	aalist=AAuniq.keys()
	aalist.sort()
		
	print
	print '#INPUT PARAMETERS'
	print 'bacground:\t',CDSfile
	print	'input file:\t',
	if LOCI!={}:
		print listfile
	else:
		print '-'
	print 'USE ONLY FIRST N CODONS:\t',
	if CUT>0:
		print CUT
	else:
		print '-'
	print
	print '#GENERAL STATISTICS'
	print
	print 'N selected (unique) IDs:\t',len(LOCI)
	print 'N good CDSs:\t',Nseqs
	print
	print 'BAD CDSs (length):\t',len(BADlen)
	print 'BAD CDSs (start):\t',len(BADstart)
	print 'BAD codons (X nts):\t',BADcod
	print 'TOT good codons:\t',TOTcod
	print
	print '%GC position 1:\t', GC1*100.0/TOTcod
	print '%GC position 2:\t', GC2*100.0/TOTcod
	print '%GC position 3:\t', GC3*100.0/TOTcod
	print '%GC TOTAL:\t', (GC1+GC2+GC3)*100.0/(TOTcod*3)
	print
	print '#CODON COUNT (per aa / 1000 codons):'
	for aa in aalist:
		print aa,'\t',AA.get(aa,0),'\t',
		print '%5.3f' % (AA.get(aa,0)*1000.0/TOTcod)
		for triplet in CODE.keys():
			if CODE[triplet]==aa:
				print triplet,'\t',TRIPLET.get(triplet,0),'\t',
				print '%5.3f' % (TRIPLET.get(triplet,0)*1000.0/TOTcod)
	print 
		