#!/usr/bin/python

#WHAT IS THIS PROGRAM FOR:
#this script extracts RNAs ID, sequences, and 2ndary structures 
#with selected constrained overlap to a given gene list, namely: 
#	- gene type = all|protein-coding
#	- feature type (only if protein-coding genes selected) = exon|CDS|5UTR|3UTR 
#	- overlap type = partial|full|full feature (otherwise full with gene but partial for selected feature)
#additionally, RNAs can be filtered by LENGTH and/or SCORE before computing the overlap itself 
#(see CONTROL PANEL below).
#GENE TYPE = ALL + OVERLAP TYPE = PARTIAL yields detailed output for characteristics of overlapping RNA/genes 
#(descriptors from original files, overlap coordinates, overlap type, and if a gene is protein-coding also detail of overlapping features)
#the processing is designed to use the RNA structure annotation from Kern&Hupalo 2013 
#(see input files below).
#the output consists of: RNA structure original ID, overlapping gene ID,  
#overlap type (partial/full, CDS,UTR,exon), coordinates, score, and length;
#an additional file is produced which contains RNA sequence and structure 
#in a format readable by RNAplot and RNAdistance (Vienna package).

#================================================

#CONTROL PANEL

SCORE=0		#INCLUSIVE, 50th percentile == 201, 95th percentile == 415	(range: 100-1217)
MINlen=0	#INCLUSIVE, minimum meaningful hairpin size == 6 (=::)		(minimum in set: 3)
MAXlen=1000	#EXCLUSIVE, change together with MINlen to create bins		(maximum in set: 730)

models=	'TAIR10_representative_gene_models.txt'
gff3=		'TAIR10_GFF3_genes.gff'
RNAannFile=	'UCSC_RNAstr_annotation.txt'
RNAseqFile=	'UCSC_RNAstr_sequence.txt'

#================================================

#load target list (if any provided)

import sys

if len(sys.argv)!=5 or int(sys.argv[2]) not in [0,1] or int(sys.argv[3]) not in [0,1,2,3,4] or int(sys.argv[4]) not in [0,1]:
	print
	print 'usage: <target list (AGI loci)> <GENEs 0-1> <OVERLAP_type 0-4> <OV.coverage 0-1>'
	print 
	print '#OPTIONS VALUES:'
	print 'GENEs:\t0 == ALL gene types, 1 == protein-coding only'
	print 'OVERLAP_type:\t0 == (with) gene, 1 == (with) exon, 2 == (with) CDS, 3 == with 5UTR, 4 == with 3UTR'
	print 'OV.coverage:\t0 == partial, 1 == full'
	print '(see script itself for RNA filtering options)'
	print '-> SET 0,0,0 TO GET DETAILED OVERLAP FIGURES FOR YOUR GENE LIST! (together with RNA structures output)'
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
	strfile=out+'.str.out'
	
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
				elif (filterG==0 and filterF==0 and filterOV==0) or (filterF==1 and feat=='exon') or (filterF==2 and feat=='CDS') or (filterF==3 and feat=='five_prime_UTR') or (filterF==4 and feat=='three_prime_UTR'):
					start,stop=int(x[3]),int(x[4])
					Gfeat[model]=Gfeat.get(model,[])+[(start,stop,feat)]
					if filterG==0 and filterF==0 and filterOV==0:
						if feat=='five_prime_UTR':
							UTR5gene[model]=1
						elif feat=='three_prime_UTR':
							UTR3gene[model]=1
							
	print '#GFF3 loaded'

#================================================

#check overlap between selected models and RNA structures

	OVERLAP={}								#gene coord. : RNA coord. (overlap according to input criteria)
	FULL,EXON,CDS,UTR5,UTR3={},{},{},{},{}	#RNA coord. : overlap with featureX (0 if partial, 1 if full) - used only if filterG,filterF,filterOV==0,0,0
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
						if filterG==0 and filterF==0 and filterOV==0:
							if (start>=startG and stop<=stopG):
								FULL[(model,R)]=1
							if Gfeat.get(model,0)!=0: #only if protein-coding gene
								for feature in Gfeat[model]:
									startF,stopF,feat=feature[0],feature[1],feature[2]
									#full overlap overwrites partial overlap with same type of feature
									if start>=startF and stop<=stopF: 
										if feat=='exon':
											EXON[(model,R)]=1
										elif feat=='CDS':
											CDS[(model,R)]=1
										elif feat=='five_prime_UTR':
											UTR5[(model,R)]=1
										else:
											UTR3[(model,R)]=1
									elif (start<=startF and stop>=startF) or (start>=startF and start<=stopF): 
								   #counted as partial only if no full overlap present for the same type of feature
										if feat=='exon' and EXON.get((model,R),-1)==-1:
											EXON[(model,R)]=0
										elif feat=='CDS' and CDS.get((model,R),-1)==-1:
											CDS[(model,R)]=0
										elif feat=='five_prime_UTR' and UTR5.get((model,R),-1)==-1:
											UTR5[(model,R)]=0
										elif UTR3.get((model,R),-1)==-1:
											UTR3[(model,R)]=0
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
	
#================================================

#load RNA sequence and structure for overlapping RNAs
	
	RNAnames={}
	for R in RNAcount.keys():
		for name in RNAcoord[R]:
			RNAnames[name]=1
	RNAstr={}
	for line in open(RNAannFile,'rU').readlines()[1:]:
		x=line.split()
		Rname= x[4]
		Rstr= x[8]
		if RNAnames.get(Rname,0)!=0:
			RNAstr[Rname]=Rstr
	RNAseq={}
	for line in open(RNAseqFile,'rU').readlines():
		x=line.split()
		if len(x)!=0:
			if x[0][0]=='>':
				READ=False
				Rname='_'.join(x[0].split('_')[2:])
				if RNAnames.get(Rname,0)!=0:
					READ=True
			elif READ:
				RNAseq[Rname]= RNAseq.get(Rname,'')+x[0]
	
	print '#RNA structure annotation extracted'
						
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
	print '#RNA SEQ.:\t',RNAseqFile
	print
	print '################################################'
	print
	print '#FILTERING OPTIONS (gene features/overlap)'
	print 
	print '#GENEs:\t', filterG 			#ALL|protein-coding only
	print '#OVERLAP type:\t',filterF	#gene feature (if protein-coding)
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
	print '#N RNA coordinates overlapping genes according to criteria:\t',len(RNAcount.keys())
	print '#(unique RNA coordinates passing filtering / total RNAs):\t',Nrna,'/',Rnames
	print 
	f=open(outfile,'w')
	s=open(strfile,'w')
	f.write('#gene_list: '+infile)
	f.write('#RNAstr.filter: S|minL|maxL : '+str(SCORE)+'|'+str(MINlen)+'|'+str(MAXlen)+'\n')
	f.write('#OVERLAPfilter: gene|feat|cover. : '+str(filterG)+'|'+str(filterF)+'|'+str(filterOV)+'\n')
	f.write('#Noverlap_genes/Nmapped_genes/Nlist_genes: '+str(len(OVsort))+'/'+str(Nloci)+'/'+str(Nlist)+'\n')
	f.write('#NoverlapRNAcoords/RNAunique_coords/TOT_RNAs: '+str(len(RNAcount.keys()))+'/'+str(Nrna)+'/'+str(Rnames)+'\n')
	f.write('#\n')
	if not (filterG==0 and filterF==0 and filterOV==0):
		f.write('#LOCUS\ttype\tchr\tstart\tstop\tstrand/S\tL\tTOT_RNAcoords\n')
		for gene in OVsort: #sorted list of gene coordinates with overlapping RNAs
			model=GENE[gene]
			f.write(model+'\t'+Gtype[model]+'\t'+gene[0]+'\t'+str(gene[1])+'\t'+str(gene[2])+'\t'+gene[3]+'\t'+str(gene[2]-gene[1]+1)+'\t'+str(len(OVERLAP[gene]))+'\n')
			for R in OVERLAP[gene]:
				for name in RNAcoord[R]:
					f.write(name+'\tRNA\t'+R[0]+'\t'+str(R[1])+'\t'+str(R[2])+'\t'+str(RNAfeat[name][0])+'\t'+str(RNAfeat[name][1])+'\t-\n')
	else:
		f.write('#LOCUS\ttype\tchr\tstart\tstop\tstrand/S\tL\tTOT_RNAcoords\tGENEp\tGENEf\tEXONp\tEXONf\tCDSp\tCDSf\t5UTRp\t5UTRf\t3UTRp\t3UTRf\n')
		for gene in OVsort: #sorted list of gene coordinates with overlapping RNAs
			model=GENE[gene]
			f.write(model+'\t'+Gtype[model]+'\t'+gene[0]+'\t'+str(gene[1])+'\t'+str(gene[2])+'\t'+gene[3]+'\t'+str(gene[2]-gene[1]+1)+'\t'+str(len(OVERLAP[gene]))+'\t')
			Gp,Gf=0,0
			for R in OVERLAP[gene]:
				if FULL.get((model,R),-1)!=-1:
					Gf+=1
				else:
					Gp+=1
			if Gfeat.get(model,0)!=0:
				Ep,Ef,Cp,Cf,p5,f5,p3,f3=0,0,0,0,0,0,0,0
				for R in OVERLAP[gene]:
					if EXON.get((model,R),-1)!=-1:
						if EXON[(model,R)]==0:
							Ep+=1
						else:
							Ef+=1
					if CDS.get((model,R),-1)!=-1:
						if CDS[(model,R)]==0:
							Cp+=1
						else:
							Cf+=1
					if UTR5.get((model,R),-1)!=-1:
						if UTR5[(model,R)]==0:
							p5+=1
						else:
							f5+=1
					if UTR3.get((model,R),-1)!=-1:
						if UTR3[(model,R)]==0:
							p3+=1
						else:
							f3+=1
			else:
				Ep,Ef,Cp,Cf,p5,f5,p3,f3='-','-','-','-','-','-','-','-'
			f.write(str(Gp)+'\t'+str(Gf)+'\t'+str(Ep)+'\t'+str(Ef)+'\t'+str(Cp)+'\t'+str(Cf)+'\t'+str(p5)+'\t'+str(f5)+'\t'+str(p3)+'\t'+str(f3)+'\n')
			for R in OVERLAP[gene]:
				if FULL.get((model,R),0)!=0:
					Gf,Gp='1','-'
				else:
					Gf,Gp='-','1'
				if Gfeat.get(model,0)!=0:
					Ep,Ef,Cp,Cf,p5,f5,p3,f3='-','-','-','-','-','-','-','-'
					if EXON.get((model,R),-1)!=-1:
						if EXON[(model,R)]==0:
							Ep='1'
						else:
							Ef='1'
					if CDS.get((model,R),-1)!=-1:
						if CDS[(model,R)]==0:
							Cp='1'
						else:
							Cf='1'
					if UTR5.get((model,R),-1)!=-1:
						if UTR5[(model,R)]==0:
							p5='1'
						else:
							f5='1'
					if UTR3.get((model,R),-1)!=-1:
						if UTR3[(model,R)]==0:
							p3='1'
						else:
							f3+='1'
				for name in RNAcoord[R]:
					f.write(name+'\tRNA\t'+R[0]+'\t'+str(R[1])+'\t'+str(R[2])+'\t'+str(RNAfeat[name][0])+'\t'+str(RNAfeat[name][1])+'\t-\t')
					if name == RNAcoord[R][0]:
						f.write(Gp+'\t'+Gf+'\t'+Ep+'\t'+Ef+'\t'+Cp+'\t'+Cf+'\t'+p5+'\t'+f5+'\t'+p3+'\t'+f3+'\n')
					else:
						f.write('\n')
						
	for gene in OVsort:
		for R in OVERLAP[gene]:
			for name in RNAcoord[R]:
				s.write('>'+name+':L'+str(RNAfeat[name][1])+'_S'+str(RNAfeat[name][0])+'\n')
				s.write(RNAseq[name]+'\n')
				s.write(RNAstr[name]+'\n')
	f.close()
	s.close()