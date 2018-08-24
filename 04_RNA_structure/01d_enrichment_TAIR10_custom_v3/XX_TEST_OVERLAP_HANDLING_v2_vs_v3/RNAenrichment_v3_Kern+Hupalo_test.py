#!/usr/bin/python

RNA='RNA_chromosome_handling.txt' #'RNA_overlap_handling.txt'
gff3='GENE_chromosome_handling.txt' #'GENE_overlap_handling.txt'

#================================================

#load RNA coordinates

RNAcoord={}
for line in open(RNA,'rU').readlines()[1:]:
	x= line.split()
	chr,start,stop,name=x[0][3:],int(x[1]),int(x[2]),x[3]
	coord=(chr,start,stop,name)
	RNAcoord[coord]=1
	
Rsort=RNAcoord.keys()
Rsort.sort()
Nrna= len(Rsort)

print '#RNA coords. loaded'

#================================================

#load gene features

Gall,Pall,TEall,miRall,snRall,snoRall,tRall,rRall,ncRall=0,0,0,0,0,0,0,0,0
UTR5genes,UTR3genes={},{}

GENE={}				#coord. : model
Gsort=[]				#sorted list of coordinates
Gtype,Gfeat={},{}	#model : type, model: features coordinates (protein-coding only)

for line in open(gff3,'rU').readlines()[1:]:
	x=line.split()
	feat=x[1]
	if feat in ['gene','pseudogene','transposable_element_gene']:
		READ=False
		chr,start,stop,model=x[0][3:],int(x[2]),int(x[3]),x[4]
		GENE[(chr,start,stop,model)]=model
		Gtype[model]=feat
		if feat=='pseudogene':
			Pall+=1
		elif feat=='transposable_element_gene':
			TEall+=1
		else: #gene
			type= x[8].split(';')[1][5:] #protein_coding_gene,other_RNA,snoRNA,snRNA,miRNA,tRNA,rRNA
			Gtype[model]=type
			if type=='protein_coding_gene':
				READ=True
				Gall+=1
			elif type=='rRNA':
				rRall+=1
			elif type=='miRNA':
				miRall+=1
			elif type=='snRNA':
				snRall+=1
			elif type=='snoRNA':
				snoRall+=1
			elif type=='tRNA':
				tRall+=1
			else: #other_RNA
				ncRall+=1
	elif READ==True and feat in ['CDS','five_prime_UTR','three_prime_UTR','exon']:
			locus=x[8].split('=')[1].split(',')[0].split()[0] #removes \n
			if locus!=model:
				continue
			else:
				start,stop=int(x[2]),int(x[3])
				Gfeat[model]=Gfeat.get(model,[])+[(start,stop,feat)]
				if feat=='five_prime_UTR':
					UTR5genes[model]=1
				elif feat=='three_prime_UTR':
					UTR3genes[model]=1

UTR5all,UTR3all=len(UTR5genes.keys()),len(UTR3genes.keys())

Gsort=GENE.keys()
Gsort.sort()
for model in Gfeat.keys(): #tested: not necessary because features already sorted in ggf3 file
	Fsort=Gfeat[model]
	Fsort.sort()
	Gfeat[model]=Fsort
Nloci=len(Gsort)

print '#GFF3 loaded'

#================================================

#check overlap between selected models and RNA structures
#(stores data that can be printed later on)

OVERLAP={}	#gene coord. : RNA coord. (any degree of overlap)
FULL={}		#gene coord. : RNA coord. (RNA fully contained in gene limits)
SCAN=0
CHR=Rsort[0][0] #initial reference chromosome for RNAs
for R in Rsort:
	FIRST=True
	chr,start,stop,name=R[0],R[1],R[2],R[3]
	for i in range(SCAN,len(Gsort)):
		gene=Gsort[i]
		chr2,start2,stop2,name2=gene[0],gene[1],gene[2],gene[3]
		print CHR,R,gene,
		if chr!=chr2: 
			print
			if chr==CHR: #new gene chromosome
				SCAN=i #keep same gene, go on with RNAs
				break
			else: #new RNA chromosome, go on with genes
				continue
		else:
			if CHR!=chr:
				CHR=chr	#update reference RNA chromosome (new chromosome)
				SCAN=i	#start scanning from 1st gene of new chromosome
			if (start<=start2 and stop>=start2) or (start>=start2 and start<=stop2):
				if FIRST==True:
					SCAN=i
				FIRST=False
				OVERLAP[gene]=OVERLAP.get(gene,[])+[R]
				if start>=start2 and stop<=stop2:
					FULL[gene]=FULL.get(gene,[])+[R]
				print '\tX!'
			elif stop<start2:	#RNA ends before gene, meaningless to increment genes (check next RNA)
				if FIRST:
					SCAN=i
				print
				break
			else: #RNA starts after gene, so might overlap with the next (check next gene)
				print
				continue

print '#OVERLAP done (general)'
			
OVsort=OVERLAP.keys()
OVsort.sort()
Fsort=FULL.keys()
Fsort.sort()

#================================================

#extracts overlap statistics (gene counts)

G,P,TE,miR,snR,snoR,tR,rR,ncR=0,0,0,0,0,0,0,0,0				#gene categories (any degree of overlap)
Gf,Pf,TEf,miRf,snRf,snoRf,tRf,rRf,ncRf=0,0,0,0,0,0,0,0,0	#gene categories (>=1 RNA fully contained in gene limits)			
c5,c3,cCDS,cEX=0,0,0,0		#protein-coding gene features (any degree of overlap)
c5f,c3f,cCDSf,cEXf=0,0,0,0	#protein-coding gene features (>=1 RNA fully contained in CDS or EXON)
for gene in OVsort:
	model=GENE[gene]
	type= Gtype[model]
	if type!='protein_coding_gene':
		if type=='pseudogene':
			P+=1
			if FULL.get(gene,0)!=0:
				Pf+=1
		elif type=='transposable_element_gene':
			TE+=1
			if FULL.get(gene,0)!=0:
				TEf+=1
		elif type=='rRNA':
			rR+=1
			if FULL.get(gene,0)!=0:
				rRf+=1
		elif type=='miRNA':
			miR+=1
			if FULL.get(gene,0)!=0:
				miRf+=1
		elif type=='snRNA':
			snR+=1
			if FULL.get(gene,0)!=0:
				snRf+=1
		elif type=='snoRNA':
			snoR+=1
			if FULL.get(gene,0)!=0:
				snoRf+=1	
		elif type=='tRNA':
			tR+=1
			if FULL.get(gene,0)!=0:
				tRf+=1		
		else: #other_RNA
			ncR+=1
			if FULL.get(gene,0)!=0:
				ncRf+=1					
	else: #protein_coding_gene
		G+=1
		if FULL.get(gene,0)!=0:
			Gf+=1
		for R in OVERLAP[gene]:
			start,stop=R[1],R[2]
			UTR5x,UTR3x,CDSx,exonx=False,False,False,False
			UTR5xf,UTR3xf,CDSxf,exonxf=False,False,False,False
			for feat in Gfeat[model]:
				start2,stop2,what=feat[0],feat[1],feat[2]
				if (start<=start2 and stop>=start2) or (start>=start2 and start<=stop2):
					if what=='five_prime_UTR':
						UTR5x=True
						if start>=start2 and stop<=stop2:
							UTR5xf=True
					elif what=='three_prime_UTR':
						UTR3x=True
						if start>=start2 and stop<=stop2:
							UTR3xf=True
					elif what=='CDS':
						CDSx=True
						if start>=start2 and stop<=stop2:
							CDSxf=True
					else: #exon
						exonx=True
						if start>=start2 and stop<=stop2:
							exonxf=True
		if exonx:
			cEX+=1
			if exonxf:
				cEXf+=1
			if UTR5x:
				c5+=1
				if UTR5xf:
					c5f+=1
			if UTR3x:
				c3+=1
				if UTR3xf:
					c3f+=1
			if CDSx:
				cCDS+=1
				if CDSxf:
					cCDSf+=1

print '#OVERLAP done (protein-coding genes)'


#================================================

#output file (counts + p-values, if background given)

print '#GENERAL STATISTICS'
print
print '#N RNAs:\t',Nrna
print '#N genes:\t',Nloci
print 
print '#N genes'
print 'ALL_GENES\t',Gall+Pall+TEall+miRall+snRall+snoRall+tRall+rRall+ncRall
print 'all_protein-coding\t',Gall
print 'all_5UTRgene\t',UTR5all
print 'all_3UTRgene\t',UTR3all
print 'all_pseudogene\t',Pall
print 'all_TrEl\t',TEall
print 'all_miRNA\t',miRall
print 'all_snRNA\t',snRall
print 'all_snoRNA\t',snoRall
print 'all_tRNA\t',tRall
print 'all_rRNA\t',rRall
print 'all_other_non-coding\t',ncRall
print
print '#N genes with RNA structure overlap (RNA partially/FULLY contained in gene limits)'
print 'protein-coding\t',G,'\t',Gf
print 'pseudogene\t',P,'\t',Pf
print 'TrEl\t',TE,'\t',TEf
print 'miRNA\t',miR,'\t',miRf
print 'snRNA\t',snR,'\t',snRf
print 'snoRNA\t',snoR,'\t',snoRf
print 'tRNA\t',tR,'\t',tRf
print 'rRNA\t',rR,'\t',rRf
print 'other_non-coding\t',ncR,'\t',ncRf
print 
print '#N protein-coding genes with RNA structure overlapping exon feature (RNA partially/FULLY contained in feature limits)'
print '5UTR\t',c5,'\t',c5f
print '3UTR\t',c3,'\t',c3f
print 'CDS\t',cCDS,'\t',cCDSf
print 'exon\t',cEX,'\t',cEXf
print
