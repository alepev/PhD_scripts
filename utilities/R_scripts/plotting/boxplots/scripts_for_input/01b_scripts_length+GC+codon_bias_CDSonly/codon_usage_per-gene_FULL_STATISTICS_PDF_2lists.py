#!/usr/bin/python

#NOTE - modified to:
#- convert loci to uppercase
#- GC & GC3 content as %, not ratio (i.e. 2 additional informative values)

import sys
import os
import subprocess as sub
sub.PIPE= 1

MIN=1	#generate counts only for genes with N codons >= MIN (after removing START e STOP codons)
			#NOTE: for Nc computations >= 100 codons recommended

if len(sys.argv)<4 or len(sys.argv)>5:
	print
	print 'usage: <CUT (0 for no cut)> <list 1> <list 2> <background (OPTIONAL)>'
	print 'compares the two lists as independent samples (default reference sequences are ORFs/CDSs);'
	print 'extracts sequence length, GC+GC3 content, and effective number of codons (Nc, Sun et al. 2012)'
	print '(both parametric (Z-test,T-test) and non-parametric (Wilcoxon) tests are used)'
	print
	print 'NOTE:'
	print '- if CUT!=0, only the first <CUT> codons (after START) will be used for the statistics'
	print '- START and STOP codons are ignored in the statistics'
	print '- bad codons (containing a character other than A,C,G,T) are also ignored'
	print '- sequences with an invalid length (!= 3N) are trimmed at the 3-end'
	print '- sequences with an invalid START or an internal STOP codon are excluded from the analysis'
	if MIN>0:
		print '- sequences with less than',MIN,'codons (except START and STOP) are also excluded from the analysis'
		print '  (see code to modify this behavior)'
	print
	
else:

	CDSfile=		'TAIR10_cds_20101214_updated.fa'
	modelfile=	'TAIR10_representative_gene_models.txt'
	
	CUT= int(sys.argv[1])
	list1= sys.argv[2]
	list2= sys.argv[3]
	bkg='-'
	if len(sys.argv)==5:
		bkg=sys.argv[4]

	#----------------------------------------------------
	# STANDARD GENETIC CODE 
	#----------------------------------------------------
	
	#(plant chloroplast and mitochondrial genomes also use the standard genetic code)
	
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
		
	#----------------------------------------------------
	# Effective Number of Codons (Nc)
	#----------------------------------------------------
	
	#revised implementation according to Sun et al. (2012), including definition of codon families below
	
	CF={} #dictionary of Codon Families, sorted by N of codons
	CF[1]=['Met','Trp']
	CF[2]={'Phe':['TTT','TTC'],'Tyr':['TAT','TAC'],'Cys':['TGT','TGC'],'Asn':['AAT','AAC'],'Gln':['CAA','CAG'],'Asp':['GAT','GAC'],'Glu':['GAA','GAG'],'Lys': ['AAA','AAG'],'His':['CAT','CAC']}
	CF[3]={'Ile':['ATT','ATC','ATA']}
	CF[4]={'Val':['GTT','GTC','GTA','GTG'],'Thr':['ACT','ACC','ACA','ACG'],'Ala':['GCT','GCC','GCA','GCG'],'Pro':['CCT','CCC','CCA','CCG'],'Gly':['GGT','GGC','GGA','GGG']}
	#families of 6 codons are split into synonimous families of 2 and 4 codons
	CF[2]['Leu']=['TTA','TTG']
	CF[2]['Ser']=['AGT','AGC']
	CF[2]['Arg']=['AGA','AGG']
	CF[4]['Leu']=['CTT','CTC','CTA','CTG']
	CF[4]['Ser']=['TCT','TCC','TCA','TCG']
	CF[4]['Arg']=['CGT','CGC','CGA','CGG']

	def SunNc(locus,COUNT,CF): 
		Nc=len(CF.get(1,[])) #N CF with 1 codon = 2 (Met,Trp)
		for m in CF.keys():
			if m!=1: 
				k= len(CF[m].keys())
				num,den=0.0,0.0
				for cf in CF[m].keys():
					Ncf,Fcf=0.0,0.0
					for triplet in CF[m][cf]:	
						Ni= COUNT[locus].get(triplet,0)	#n codons
						Ncf+= Ni #tot codons in synonimous family 
						Fcf+= (Ni+1)**2
					num+= Ncf
					den+= Ncf*Fcf/(Ncf+m)**2
				if den!=0:
					Nc+= k*num/den
				else:
					#print locus (< 0.1% of loci seem to lack a cf, so very minor problem)
					Nc+= 0 #if a cf is not used at all (most likely this could happen for Ile, i.e. only representative of cf with 3 codons)
					#alternatively, one could assume equal usage for a missing cf by using Nc+= k*m)
		return Nc
		
	#----------------------------------------------------
	# GC & GC3 content
	#----------------------------------------------------
	
	#GC excludes START and STOP, GC3s additionally excludes 1-codon families (Met and Trp)
	
	def GCcount(locus,COUNT,CODE): 
		GCs,GC3s=0,0
		codons,codons2=0,0
		for triplet in COUNT[locus].keys():
			N=COUNT[locus][triplet]
			codons+=N
			if CODE[triplet] not in ['Met','Trp']:
					codons2+=N
			if triplet[0] in ['G','C']:
				GCs+=N
			if triplet[1] in ['G','C']:
				GCs+=N
			if triplet[2] in ['G','C']:
				GCs+=N
				if CODE[triplet] not in ['Met','Trp']:
					GC3s+=N
		return GCs*100.0/(3*codons), GC3s*100.0/codons2
	
	#----------------------------------------------------
	# TEST list1 vs list2
	#----------------------------------------------------
	
	#calculates Z-test + T-test + Wilcoxon non-parametric test p-values
	#generates boxplot and kernel density plot of the two distributions + mean and median annotation

	def test(valuedict,locidict,name,VALUE,bkg): #bkg here is a dictionary of loci in the background population
	
		tmp1=name+'.list1.tmp'
		tmp2=name+'.list2.tmp'
		pdf1=name+'.density.pdf'
		pdf2=name+'.boxplot.pdf'
		
		if bkg!='-':
			tmp3=name+'.bkg.tmp'
		else:
			tmp3=bkg
	
		l1,l2=0,0
		
 		f1=open(tmp1,'w')
 		f2=open(tmp2,'w')
 		if bkg!='-':
			f3=open(tmp3,'w')
 		for locus in valuedict.keys():
 			if locidict.get(locus,-1)!=-1:
 				if locidict[locus]<=1:
  					f1.write(str(valuedict[locus])+'\n')
  					l1+=1
  				if locidict[locus]>=1:
  					f2.write(str(valuedict[locus])+'\n')
  					l2+=1
  			if bkg!='-' and bkg.get(locus,0)!=0:
  				f3.write(str(valuedict[locus])+'\n')
 		f1.close()
 		f2.close()
 		if bkg!='-':
			f3.close()
		
 		x=sub.Popen('Rscript test_list1_vs_list2_PDF.R '+tmp1+' '+tmp2+' '+pdf1+' '+pdf2+' '+VALUE+' '+tmp3, shell=True, stdout=1)
 		out= x.communicate()[0].split() 
 		os.remove(tmp1)
		os.remove(tmp2)
		if bkg!='-':
			os.remove(tmp3)
		
		m1= float(out[0])
		m2= float(out[1])
		sd1= float(out[2])
		sd2= float(out[3])
		me1= float(out[4])
		me2= float(out[5])
		iqr1= float(out[6])
		iqr2= float(out[7])
		Zpval= float(out[8])
		Tpval= float(out[9])
		Wpval= float(out[10])
		if bkg!='-':
			popvar= float(out[11])
		else:
			popvar='-'

		print '#list1'
		print 'N:\t\t',l1
		print 'mean (sd):\t%1.2f (%1.2f)' % (m1,sd1)
		print 'median (IQR):\t%1.2f (%1.2f)' % (me1,iqr1)
		print '#list2'
		print 'N:\t\t',l2
		print 'mean (sd):\t%1.2f (%1.2f)' % (m2,sd2)
		print 'median (IQR):\t%1.2f (%1.2f)' % (me2,iqr2)
		print
		print 'Z-test p-value (N>30):\t%1.2e' % Zpval
		if bkg!='-':
			print '(background pop. variance provided)'
		else:
			print '(background pop. variance not available)'
		print 'T-test p-value (N<30):\t%1.2e' % Tpval
		print 'Wilcoxon test p-value:\t%1.2e' % Wpval
		print
		print '#check density plot in:\t',pdf1
		print '#check boxplot in:\t',pdf2
		print
		
	#----------------------------------------------------
	# GENERAL FUNCTIONS
	#----------------------------------------------------
	
	#generate list of keys for which dict[key]==value
	def member(dict,value):
		return filter(lambda x:dict[x]==value,dict)
		
	#generate list of keys for which dict[key] <= or >= value
	def membercfr(dict,value,cfr):
		if cfr=='<=':
			return filter(lambda x:dict[x]<=value,dict)
		elif cfr=='>=':
			return filter(lambda x:dict[x]>=value,dict)
	
	#split sequence/string into a list of words with specified size
	def splitbins(seq,wordsize):
		if len(seq)%wordsize==0:
			return [seq[i:i+wordsize] for i in range(0, len(seq), wordsize)]
		else:
			return [seq[i:i+wordsize] for i in range(0, len(seq), wordsize)]+[seq[len(seq)-len(seq)%wordsize:]]
    	
    	#in practice the latter is never used if sequence are pre-filtered for this
				
	#####################################################
	# PROGRAM START
	#####################################################

	#load representative gene models

	MOD={}
	for line in open(modelfile,'rU').readlines():
			x=line.split()
			if len(x)!=0 and x[0][:2]=='AT':
				locus=x[0].split('.')
				MOD[locus[0]]=locus[1]
				
	#load list1
	
	NOTFOUND={}
	LOCI={} #locus in list1: 0, list1+2: 1, list2: 2
	for line in open(list1,'rU').readlines():
		x=line.split()[0].upper()
		if len(x)>0 and x[0]!='#':
			if x[0]=='>':
				x=x[1:]
			if x[:2]=='AT':
				locus=x.split('|')[0].split('/')[0] #remove descriptors, add here other separators if needed
				if len(locus.split('.'))>1: #gene model provided
					LOCI[locus]=0
				elif MOD.get(locus,0)!=0: #representative gene model automatically assigned
					LOCI[locus+'.'+MOD[locus]]=0
				else:
					NOTFOUND[locus]=0
				
	#load list2
	
	for line in open(list2,'rU').readlines():
		x=line.split()[0].upper()
		if len(x)>0 and x[0]!='#':
			if x[0]=='>':
				x=x[1:]
			if x[:2]=='AT':
				locus=x.split('|')[0].split('/')[0] #remove descriptors, add here other separators if needed
				if len(locus.split('.'))>1: #gene model provided
					if LOCI.get(locus,-1)==0:
						LOCI[locus]=1 #also in list1
					elif LOCI.get(locus,-1)!=1:
						LOCI[locus]=2
				elif MOD.get(locus,0)!=0: #representative gene model automatically assigned
					locus=locus+'.'+MOD[locus]
					if LOCI.get(locus,-1)==0:
						LOCI[locus]=1 #also in list1
					elif LOCI.get(locus,-1)!=1:
						LOCI[locus]=2
				else:
					if NOTFOUND.get(locus,-1)==0:
						NOTFOUND[locus]=1 #also in list1
					else:
						NOTFOUND[locus]=2
						
	#load background population (if provided)
	
	BKG='-'
	if bkg!='-':
		BKG={} 
		for line in open(bkg,'rU').readlines():
			x=line.split()[0].upper()
			if len(x)>0 and x[0]!='#':
				if x[0]=='>':
					x=x[1:]
				if x[:2]=='AT':
					locus=x.split('|')[0].split('/')[0] #remove descriptors, add here other separators if needed
					if len(locus.split('.'))>1: #gene model provided
						BKG[locus]=1
					elif MOD.get(locus,0)!=0: #representative gene model automatically assigned
						BKG[locus+'.'+MOD[locus]]=1
		
	#scan CDS file for codon counts
	
	STOPcodons=member(CODE,'STOP') #used for check of internal stop codons

	COUNT={}
	LEN={}			#(used only if CUT==0)
	BADlen3={}		#included with warnings (trimmed)
	BADstopEND={}	#included with warnings
	BADstart={}		#excluded
	BADlenMIN={}	#excluded
	BADstop={}		#excluded (i.e. STOP codon within CDS, before the end)
	BADcod={}		#included with warning (exclude individual codon)
	TOTcod1,TOTcod2=0,0
	TOTbad1,TOTbad2=0,0
	Nseqs1,Nseqs2=0,0
	locus=''
	seq=''
	READ=False
	for line in open(CDSfile,'rU').readlines():
			if line[0]=='>':
				if locus!='' and READ: #if new ID, check previous sequence
					#sequence checks
					OK=True
					if seq[-1]=='*':
						seq=seq[:-1]
					if len(seq)%3!=0 and LOCI.get(locus,-1)!=-1:						#used with warnings
						BADlen3[locus]=LOCI[locus]
						seq=seq[:len(seq)-len(seq)%3]
					if CODE.get(seq[-3:],0)!='STOP' and LOCI.get(locus,-1)!=-1:	#used with warnings
						BADstopEND[locus]=LOCI[locus]
					else:
						seq=seq[:-3]
					if CODE.get(seq[:3])!='Met' and LOCI.get(locus,-1)!=-1:		#excluded
						BADstart[locus]=LOCI[locus]
						OK=False
					else:
						seq=seq[3:]
					if len(seq)<MIN*3 and LOCI.get(locus,-1)!=-1:					#excluded
						BADlenMIN[locus]=LOCI[locus]
						OK=False
					if any(map(lambda x: x in STOPcodons, splitbins(seq,3))) and LOCI.get(locus,-1)!=-1: #excluded
						BADstop[locus]=LOCI[locus]
						OK=False
					#codon counts for sequences passing checks
					if OK:
						COUNT[locus]={}
						if CUT==0:
							MAX=len(seq)
							LEN[locus]=MAX
						else:
							MAX=min(len(seq),CUT*3)
						for i in range(0,MAX,3):
							triplet=seq[i:i+3]
							if CODE.get(triplet,0)!=0:
								COUNT[locus][triplet]= COUNT[locus].get(triplet,0)+1
								if LOCI.get(locus,-1)!=-1:
									if LOCI[locus]<=1:
										TOTcod1+=1
									if LOCI[locus]>=1:
										TOTcod2+=1
							elif LOCI.get(locus,-1)!=-1:
								BADcod[locus]= LOCI[locus] #used with warnings
								if LOCI[locus]<=1:
								 TOTbad1+=1
								if LOCI[locus]>=1:
									TOTbad2+=1
											
				locus=line.split()[0][1:] #save new sequence ID
				seq=''
				if LOCI.get(locus,-1)!=-1:
					if LOCI[locus]<=1:
						Nseqs1+=1
					if LOCI[locus]>=2:
						Nseqs2+=1
					READ=True
				elif bkg!='-' and BKG.get(locus,0)!=0:
					READ=True
				else:
					READ=False			
					
			elif READ: #if sequence and corresponding ID was OK, read it
				if seq=='':
					seq=line.split()[0].upper()
				else:
					seq+=line.split()[0].upper() #add to existing
						
	#last sequence in the file (see above for description)
	if locus!='' and READ: #if new ID, check previous sequence
		#sequence checks
		OK=True
		if seq[-1]=='*':
			seq=seq[:-1]
		if len(seq)%3!=0 and LOCI.get(locus,-1)!=-1:						#used with warnings
			BADlen3[locus]=LOCI[locus]
			seq=seq[:len(seq)-len(seq)%3]
		if CODE.get(seq[-3:],0)!='STOP' and LOCI.get(locus,-1)!=-1:	#used with warnings
			BADstopEND[locus]=LOCI[locus]
		else:
			seq=seq[:-3]
		if CODE.get(seq[:3])!='Met' and LOCI.get(locus,-1)!=-1:		#excluded
			BADstart[locus]=LOCI[locus]
			OK=False
		else:
			seq=seq[3:]
		if len(seq)<MIN*3 and LOCI.get(locus,-1)!=-1:					#excluded
			BADlenMIN[locus]=LOCI[locus]
			OK=False
		if any(map(lambda x: x in STOPcodons, splitbins(seq,3))) and LOCI.get(locus,-1)!=-1: #excluded
			BADstop[locus]=LOCI[locus]
			OK=False
		#codon counts for sequences passing checks
		if OK:
			COUNT[locus]={}
			if CUT==0:
				MAX=len(seq)
				LEN[locus]=MAX
			else:
				MAX=min(len(seq),CUT*3)
			for i in range(0,MAX,3):
				triplet=seq[i:i+3]
				if CODE.get(triplet,0)!=0:
					COUNT[locus][triplet]= COUNT[locus].get(triplet,0)+1
					if LOCI.get(locus,-1)!=-1:
						if LOCI[locus]<=1:
							TOTcod1+=1
						if LOCI[locus]>=1:
							TOTcod2+=1
				elif LOCI.get(locus,-1)!=-1:
					BADcod[locus]= LOCI[locus] #used with warnings
					if LOCI[locus]<=1:
					 TOTbad1+=1
					if LOCI[locus]>=1:
						TOTbad2+=1
									
	#per-gene Nc, GC, GC3 statistics
	
	Nc,GC,GC3={},{},{}
	for locus in COUNT.keys():
		Nc[locus]=SunNc(locus,COUNT,CF)
		gcount=GCcount(locus,COUNT,CODE)
		GC[locus]=gcount[0]
		GC3[locus]=gcount[1]
							
	#####################################################
	# OUTPUT
	#####################################################
				
	#sequence statistics
		
	print
	print '#INPUT PARAMETERS'
	print 'reference sequences:\t',CDSfile
	print 'list1:\t',list1
	print	'list2:\t',list2
	print '- background population (Z-test variance):\t',bkg
	print '- MIN required nts:\t',
	if MIN>0:
		print MIN
	else:
		print '-'
	print '- MAX codons after START:\t',
	if CUT>0:
		print CUT
	else:
		print '-'
	print '(n.b. codon counts exclude START and STOP)'
	print
	print '#SEQUENCE STATISTICS (LIST1)'
	print 'loci with valid model:\t',len(membercfr(LOCI,1,'<='))
	print 'of which with valid sequence:\t',Nseqs1
	print 'of which passing filtering:\t',len(filter(lambda x:COUNT.get(x,0)!=0,membercfr(LOCI,1,'<=')))
	print '#loci with warnings: (NOTE: overlap is possible)'
	print '- NOT FOUND: loci without model:\t',len(membercfr(NOTFOUND,1,'<='))
	print '- included: length not 3x:\t',len(membercfr(BADlen3,1,'<='))
	print '- included: missing STOP:\t',len(membercfr(BADstopEND,1,'<='))
	print '- included: bad codons:\t',len(membercfr(BADcod,1,'<='))
	print '- EXCLUDED: bad START:\t',len(membercfr(BADstart,1,'<='))
	print '- EXCLUDED: internal STOP:\t',len(membercfr(BADstop,1,'<='))
	print '- EXCLUDED: below MIN length:\t',len(membercfr(BADlenMIN,1,'<='))
	print '#codon counts:'
	print 'N valid codons:\t',TOTcod1
	print 'N bad codons - EXCLUDED:\t',TOTbad1
	print
	print '#SEQUENCE STATISTICS (LIST2)'
	print 'loci with valid model:\t',len(membercfr(LOCI,1,'>='))
	print 'of which with valid sequence:\t',Nseqs2
	print 'of which passing filtering:\t',len(filter(lambda x:COUNT.get(x,0)!=0,membercfr(LOCI,1,'>=')))
	print '#loci with warnings: (NOTE: overlap is possible)'
	print '- NOT FOUND: loci without model:\t',len(membercfr(NOTFOUND,1,'>='))
	print '- included: length not 3x:\t',len(membercfr(BADlen3,1,'>='))
	print '- included: missing STOP:\t',len(membercfr(BADstopEND,1,'>='))
	print '- included: bad codons:\t',len(membercfr(BADcod,1,'>='))
	print '- EXCLUDED: bad START:\t',len(membercfr(BADstart,1,'>='))
	print '- EXCLUDED: internal STOP:\t',len(membercfr(BADstop,1,'>='))
	print '- EXCLUDED: below MIN length:\t',len(membercfr(BADlenMIN,1,'>='))
	print '#codon counts:'
	print 'N valid codons:\t',TOTcod2
	print 'N bad codons - EXCLUDED:\t',TOTbad2
	
	#LIST1 vs LIST2

	name= '.'.join(list1.split('.')[:-1])+'_vs_'+'.'.join(list2.split('.')[:-1])
	if CUT>0:
		name+='.CUT'+str(CUT)
	print
	print '#LIST1 vs LIST2'
	print
	print '#effective number of codons (Nc)'
	test(Nc,LOCI,name+'.Nc','"Nc index"',BKG) 
	print
	print '#GC content'
	test(GC,LOCI,name+'.GC','"GC%"',BKG)
	print
	print '#GC3 content'
	test(GC3,LOCI,name+'.GC3','"GC3%"',BKG)
	if CUT==0:
		print
		print '#sequence length'
		test(LEN,LOCI,name+'.LEN','"sequence length (nt)"',BKG)
	print
		