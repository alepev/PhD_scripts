#!/usr/bin/python

#WHAT IS THIS PROGRAM FOR:
#this script extracts a list of SCORE and LENGTH pairs for use with a separate R script (to plot density distributions)
#for each RNA satisfying the overlap conditions with a given gene list, namely: 
#	- gene type = all|protein-coding
#	- feature type (only if protein-coding genes selected) = exon|CDS|5UTR|3UTR 
#	- overlap type = partial|full|full feature (otherwise full with gene but partial for selected feature)
#additionally, RNAs can be filtered by LENGTH and/or SCORE before computing the overlap itself 
#(see CONTROL PANEL below).
#GENE TYPE = ALL + OVERLAP TYPE = PARTIAL yields detailed output for characteristics of overlapping RNA/genes 
#(descriptors from original files, overlap coordinates, overlap type, and if a gene is protein-coding also detail of overlapping features)
#the processing is designed to use the RNA structure annotation from Kern&Hupalo 2013 
#(see input files below).

#================================================

#CONTROL PANEL

SCORE=0		#INCLUSIVE, 50th percentile == 201, 95th percentile == 415	(range: 100-1217 		| DEFAULT: 0)
MINlen=0		#INCLUSIVE, minimum meaningful hairpin size == 6 (=::)		(minimum in set: 3 	| DEFAULT: 0)
MAXlen=1000	#EXCLUSIVE, change together with MINlen to create bins		(maximum in set: 730	| DEFAULT: 1000)

models=	'TAIR10_representative_gene_models.txt'
gff3=		'TAIR10_GFF3_genes.gff'
RNAannFile=	'UCSC_RNAstr_annotation.txt'

#================================================

#load target list (if any provided)

import sys

if len(sys.argv)!=5 or int(sys.argv[2]) not in [0,1] or int(sys.argv[3]) not in [0,1,2,3,4] or int(sys.argv[4]) not in [0,1]:
	print
	print 'usage: <target list (AGI loci)> <GENEs 0-1> <OVERLAP_type 0-4> <OV.coverage 0-1>'
	print 
	print '#OPTIONS VALUES:'
	print 'GENEs:\t0 == ALL gene types, 1 == protein-coding only'
	print 'OVERLAP_type:\t0 == gene, 1 == exon, 2 == CDS, 3 == 5UTR, 4 == 3UTR'
	print 'OV.coverage:\t0 == partial, 1 == full'
	print '(see script itself for RNA filtering options)'
	print
	print '#NOTE:'
	print 'if GENEs == 0, then OV.type is automatically set to 0 (features overlap computed only if protein-coding genes are selected)'
	print
else:
	print
	print '#PROGRESS REPORT'
	print
	infile= sys.argv[1] 
	filterG= int(sys.argv[2])
	if filterG==0:
		filterF=0
	else:
		filterF= int(sys.argv[3])
	filterOV= int(sys.argv[4])
	list={}
	for line in open(infile,'rU').readlines():
		x= line.split()
		if len(x)!=0 and x[0][0]!='#':
			locus= x[0]
			list[locus]=1
			
	out='.'.join(infile.split('.')[:-1])+'.'+str(filterG)+str(filterF)+str(filterOV)
	if MINlen>0:
		out+='.L'+str(MINlen)
	if MAXlen<1000:
		out+='-'+str(MAXlen)
	if SCORE>0:
		out+='S'+str(SCORE)
	outfile=out+'.out'
	
	Nlist=len(list.keys())
	print '#gene list (AGI loci) loaded'

#================================================

#load representative gene models

	MOD={}
	for line in open(models,'rU').readlines()[3:]:
		x=line.split()
		if len(x)>0:
			locus= x[0].split('.')
			if list.get(locus[0],0)!=0:
				MOD[locus[0]]=locus[1] #locus AGI code + representative model
			
	print '#representative models loaded'
	Nloci= len(MOD.keys()) 
	if Nloci<Nlist:
		print
		print '***WARNING:',Nlist-Nloci,'provided loci with MISSING MODEL (excluded from further analyses):***'
		for locus in list.keys():
			if MOD.get(locus,0)==0:
				print locus
		print
						
#================================================

#load RNA structures coordinates

	RNAcoord={}
	RNAfeat={}
	for line in open(RNAannFile,'rU').readlines()[1:]:
		x= line.split()
		chr,start,stop,name,score,L=x[1][3:],int(x[2])+1,int(x[3])+1,x[4],int(x[5]),int(x[7])
		#NOTES: 
		# - no strand info necessary since all on +
		# - coordinates starting at 0, +1 to allow comparison to gff annotation (coord.s starting at 1)
		if score>=SCORE and L>=MINlen and L<MAXlen:
			coord=(chr,start,stop)
			RNAcoord[coord]= RNAcoord.get(coord,[])+[name]
			RNAfeat[name]=(score,L)
		
	Rsort=RNAcoord.keys()
	Rsort.sort()
	Nrna= len(Rsort)
	Rnames=len(RNAfeat.keys())
	
	print '#RNA loaded'

#================================================

#load gene features

	GENE={}						#coord. : model
	Gsort=[]						#sorted list of gene coordinates
	Gtype,Gfeat={},{}			#model : type, model : features coordinates (protein-coding only)
	UTR5gene,UTR3gene={},{}	#model : has UTR (1 if yes) - used only if filterG,filterF,filterOV==0,0,0 to keep track of protein-coding genes with UTRs
	
	for line in open(gff3,'rU').readlines()[1:]:
		x=line.split()
		feat=x[2]
		if feat in ['gene','pseudogene','transposable_element_gene']:
			READ=False
			locus=x[8].split(';')[0][3:]
			if MOD.get(locus,0)==0:
				continue
			else:
				type= x[8].split(';')[1][5:] #protein_coding_gene,other_RNA,snoRNA,snRNA,miRNA,tRNA,rRNA
				if filterG==0 or type=='protein_coding_gene': #filterG==1
					model=locus+'.'+MOD[locus] #right model
					chr,start,stop,strand=x[0][3:],int(x[3]),int(x[4]),x[6]
					GENE[(chr,start,stop,strand)]=model
					Gsort+= [(chr,start,stop,strand)]
					if feat=='gene':
						Gtype[model]=type
						if type=='protein_coding_gene':
							READ=True
					else:
						Gtype[model]=feat
		elif READ==True and feat in ['exon','CDS','five_prime_UTR','three_prime_UTR']:
				locus=x[8].split('=')[1].split(',')[0].split()[0] #removes \n
				if locus!=model:
					continue
				#information on features is stored also when input== 0,0,0 to provide detailed output
				elif (filterF==1 and feat=='exon') or (filterF==2 and feat=='CDS') or (filterF==3 and feat=='five_prime_UTR') or (filterF==4 and feat=='three_prime_UTR'):
					start,stop=int(x[3]),int(x[4])
					Gfeat[model]=Gfeat.get(model,[])+[(start,stop,feat)]
							
	print '#GFF3 loaded'

#================================================

#check overlap between selected models and RNA structures

	OVERLAP={}										#gene coord. : RNA coord. (overlap according to input criteria)
	RNAcount={}
	SCAN=0
	for R in Rsort: #RNA coordinates, might correspond to multiple RNA structures
		FIRST=True
		chr,start,stop=R[0],R[1],R[2]
		for i in range(SCAN,len(Gsort)):
			gene=Gsort[i]
			model=GENE[gene]
			chrG,startG,stopG=gene[0],gene[1],gene[2]
			if chr!=chrG:
				SCAN=i
				break
			elif (start<=startG and stop>=startG) or (start>=startG and start<stopG):
				if FIRST==True:
					SCAN=i
				FIRST=False
				if filterF==0: #whole gene
					if filterOV==0 or (start>=startG and stop<=stopG): #filterOV==1
						#either full overlap, or partial overlap allowed (gene)
						if R not in OVERLAP.get(gene,[]):
							OVERLAP[gene]=OVERLAP.get(gene,[])+[R]	
							RNAcount[R]=1
				elif Gfeat.get(model,0)!=0: #filter for specific protein-coding gene features
					for feature in Gfeat[model]:
						startF,stopF,feat=feature[0],feature[1],feature[2]
						if (start>=startF and stop<=stopF) or (filterOV==0 and ((start<=startF and stop>=startF) or (start>=startF and start<=stopF))):
						#either full overlap, or partial overlap allowed (feature)
							if R not in OVERLAP.get(gene,[]):
								OVERLAP[gene]=OVERLAP.get(gene,[])+[R]
								RNAcount[R]=1
			elif stop<startG: #RNA ends before gene, meaningless to continue
				break
			else: #RNA starts after gene, so might overlap with the next
				continue
	
	print '#OVERLAP done'
				
	OVsort=OVERLAP.keys()
	OVsort.sort()
	Rstr=0
	for R in RNAcount.keys():
		Rstr+=len(RNAcoord[R])
						
#================================================

#output summary + RNA structure file (RNAfold format)

	print 
	print '################################################'
	print
	print '#INPUT FILES'
	print 
	print '#INPUT LIST:\t',infile
	print
	print '#TAIR MODELS:\t',models
	print '#GFF3 ANNOT.:\t',gff3
	print '#RNA ANNOT.:\t',RNAannFile
	print
	print '################################################'
	print
	print '#FILTERING OPTIONS (gene features/overlap)'
	print 
	print '#GENEs:\t', filterG 			#ALL|protein-coding only
	print '#OVERLAP type:\t',filterF		#gene feature (if protein-coding)
	print '#OVERLAP cov.:\t',filterOV	#full/partial
	print 
	print '################################################'
	print
	print '#FILTERING OPTIONS (RNA structure)'
	print 
	print '#min SCORE:\t',SCORE
	print '#MIN length:\t',MINlen
	print '#MAX length:\t',MAXlen
	print 
	print '################################################'
	print
	print '#OUTPUT STATISTICS'
	print
	print '#N genes overlapping RNA coordinates according to criteria:\t',len(OVsort)
	print '#(N genes with mapped model / original genes):\t',Nloci,'/',Nlist
	print '#N RNA structures/unique coordinates overlapping genes according to criteria:\t',Rstr,'/',len(RNAcount.keys())
	print '#(unique RNA coordinates passing filtering / total coordinates):\t',Nrna,'/',Rnames
	print 
	f=open(outfile,'w')
	f.write('LIST\tSCORE\tLENGTH\tNAME\n')
	LS= infile.split('/')[-1].split('.')[0]+'.'+str(filterG)+str(filterF)+str(filterOV)
	if (SCORE!=0 or MINlen!=0 or MAXlen!=1000):
		LS+='.'
		if SCORE!=0:
			LS+='S'+str(SCORE)
		if MINlen!=0:
			LS+='L'+str(MINlen)
		if MAXlen!=1000:
			LS+='-'+str(MAXlen)
	LS+= ' (N='+str(Rstr)+')'
	for gene in OVsort: #sorted list of gene coordinates with overlapping RNAs
		for R in OVERLAP[gene]:
				for name in RNAcoord[R]: 
					f.write(LS+'\t'+str(RNAfeat[name][0])+'\t'+str(RNAfeat[name][1])+'\t'+name+'\n')
	f.close()