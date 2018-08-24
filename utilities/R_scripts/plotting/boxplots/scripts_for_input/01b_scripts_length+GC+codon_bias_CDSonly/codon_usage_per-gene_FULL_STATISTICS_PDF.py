#!/usr/bin/python

#NOTE - modified to:
#- convert loci to uppercase
#- NOT remove temporary lists when single file (for re-openining and plotting in R)
#- GC & GC3 content as %, not ratio (i.e. 2 additional informative values)

import sys
import os
import subprocess as sub
sub.PIPE= 1

MIN=100	#generate counts only for genes with N codons >= MIN (after removing START e STOP codons)
			#NOTE: for Nc computations >= 100 codons recommended

if len(sys.argv)<3 or len(sys.argv)>4:
	print
	print 'usage: <CUT (0 for no cut)> <list 1> <list 2 (OPTIONAL)>'
	print 'extracts sequence length, GC+GC3 content, and effective number of codons (Nc, Sun et al. 2012)'
	print 'for each gene after mapping it to corresponding CDS;'
	print '- if only list1, tests for normality of the distributions (qqplot, histogram, and normality test)'
	print '- if both list1 and list2, assumes list1 as population and tests for significantly different mean/median of list2'
	print '  (uses both parametric (Z-test,T-test) and non-parametric (Wilcoxon) tests,'
	print '   see output NORMALITY CHECK plots to decide what test p-value result to look at)'
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
	list1= sys.argv[2] #background, ex. genes represented on the microarray
	list2='-'
	if len(sys.argv)==4:
		list2=sys.argv[3]

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
	# NORMALITY CHECK
	#----------------------------------------------------
	
	#calculates Shapiro-Wilk p-value
	#generates Q-Q plot and histogram (100 bins) with superimposed normal distribution
	
	def normal(locidict,usekeys,name):
	
		alpha=0.05
	
		tmp=name+'.tmp'
		pdf1=name+'.QQplot.pdf'
		pdf2=name+'.histo.pdf'
		
		f=open(tmp,'w')
		for locus in usekeys:
			f.write(str(locidict[locus])+'\n')
		f.close()
		
		x=sub.Popen('Rscript test_normality_PDF.R '+tmp+' '+pdf1+' '+pdf2+' ', shell=True, stdout=1)
		out= x.communicate()[0].split() 
		#os.remove(tmp)
		
		m=  float(out[0])
		sd= float(out[1])
		pval= float(out[2]) #Shapiro-Wilk test if N<5000, Kolmogorov-Smirnov otherwise
		
		print 'mean (sd):\t%1.2f (%1.2f)' % (m,sd)
		if len(usekeys)<=5000:
			print 'Shapiro-Wilk test (N<=5000): distribution is',
		else:
			print 'Kolmogorov-Smirnov test (N>5000): distribution is',
		if pval <= alpha:
			print 'NOT NORMAL',
		else:
			print 'NORMAL',
		print '(p = %1.2e , (alpha = %1.2f) )' % (pval,alpha)
		print 
		print '#check qqplot in:\t',pdf1
		print '#check histogram in:\t',pdf2
	
	#----------------------------------------------------
	# TEST list vs background (when list1 + list2 provided)
	#----------------------------------------------------
	
	#calculates Z-test + T-test + Wilcoxon non-parametric test p-values
	#generates boxplot and kernel density plot of the two distributions + mean and median annotation

	def test(locidict,list2,name,VALUE):
	
		tmp1=name+'.bkg.tmp'
		tmp2=name+'.list.tmp'
		pdf1=name+'.density.pdf'
		pdf2=name+'.boxplot.pdf'
		
 		f=open(tmp1,'w')
 		for locus in locidict.keys():
  			f.write(str(locidict[locus])+'\n')
 		f.close()
 		
 		f=open(tmp2,'w')
 		for locus in list2:
 			f.write(str(locidict[locus])+'\n')
 		f.close()
 		
 		x=sub.Popen('Rscript test_list_vs_bkg_PDF.R '+tmp1+' '+tmp2+' '+pdf1+' '+pdf2+' '+VALUE+' ', shell=True, stdout=1)
 		out= x.communicate()[0].split() 
 		os.remove(tmp1)
		os.remove(tmp2)
		
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

		print '#background (list1)'
		print 'N:\t\t',len(locidict.keys())
		print 'mean (sd):\t%1.2f (%1.2f)' % (m1,sd1)
		print 'median (IQR):\t%1.2f (%1.2f)' % (me1,iqr1)
		print '#sample (list2)'
		print 'N:\t\t',len(list2)
		print 'mean (sd):\t%1.2f (%1.2f)' % (m2,sd2)
		print 'median (IQR):\t%1.2f (%1.2f)' % (me2,iqr2)
		print
		print 'Z-test p-value (N>30):\t%1.2e' % Zpval
		print 'T-test p-value (N<30):\t%1.2e' % Tpval
		print 'Wilcoxon test p-value:\t%1.2e' % Wpval
		print
		print '#check density plot output in:\t',pdf1
		print '#check boxplot in:\t\t',pdf2

	#----------------------------------------------------
	# GENERAL FUNCTIONS
	#----------------------------------------------------
	
	#generate list of keys for which dict[key]==value
	def member(dict,value):
		return filter(lambda x:dict[x]==value,dict)
	
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
				
	#load list1 (if list2 provided, this is the background)
	
	LOCI={}
	NOTFOUND={}
	for line in open(list1,'rU').readlines():
			x=line.split()[0].upper()
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
				
	#load list2 (if provided, this is the sample to test against the background)
	
	#(if locus in list1 only LOCI[locus]==1, if in list1+list2 LOCI[locus]==2)

	if len(sys.argv)==4:
		for line in open(list2,'rU').readlines():
			x=line.split()[0].upper()
			if len(x)>0 and x[0]!='#':
				if x[0]=='>':
					x=x[1:]
				if x[:2]=='AT':
					locus=x.split('|')[0].split('/')[0]
					if len(locus.split('.'))>1 and LOCI.get(locus,0)!=0:
							LOCI[locus]=2
					elif MOD.get(locus,0)!=0 and LOCI.get(locus+'.'+MOD[locus],0)!=0:
						LOCI[locus+'.'+MOD[locus]]=2
					else:
						NOTFOUND[locus]=2
		
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
	TOTcod,TOTcod2=0,0
	TOTbad,TOTbad2=0,0
	Nseqs,Nseqs2=0,0
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
					if len(seq)%3!=0:						#used with warnings
						BADlen3[locus]=LOCI[locus]
						seq=seq[:len(seq)-len(seq)%3]
					if CODE.get(seq[-3:],0)!='STOP':	#used with warnings
						BADstopEND[locus]=LOCI[locus]
					else:
						seq=seq[:-3]
					if CODE.get(seq[:3])!='Met':		#excluded
						BADstart[locus]=LOCI[locus]
						OK=False
					else:
						seq=seq[3:]
					if len(seq)<MIN*3:					#excluded
						BADlenMIN[locus]=LOCI[locus]
						OK=False
					if any(map(lambda x: x in STOPcodons, splitbins(seq,3))):
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
								TOTcod+=1
								if LOCI[locus]==2:
									TOTcod2+=1
							else:
								BADcod[locus]= LOCI[locus] #used with warnings
								TOTbad+=1
								if LOCI[locus]==2:
									TOTbad2+=1
											
				locus=line.split()[0][1:] #save new sequence ID
				seq=''
				if LOCI.get(locus,0)!=0:
					Nseqs+=1
					if LOCI.get(locus,0)==2:
						Nseqs2+=1
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
		OK=True
		if seq[-1]=='*':
			seq=seq[:-1]
		if len(seq)%3!=0:						#used with warnings
			BADlen3[locus]=LOCI[locus]
			seq=seq[:len(seq)-len(seq)%3]
		if CODE.get(seq[-3:],0)!='STOP':	#used with warnings
			BADstopEND[locus]=LOCI[locus]
		else:
			seq=seq[:-3]
		if CODE.get(seq[:3])!='Met':		#excluded
			BADstart[locus]=LOCI[locus]
			OK=False
		else:
			seq=seq[3:]
		if len(seq)<MIN*3:					#excluded
			BADlenMIN[locus]=LOCI[locus]
			OK=False
		if any(map(lambda x: x in STOPcodons, splitbins(seq,3))):
			BADstop[locus]=LOCI[locus]
			OK=False
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
					TOTcod+=1
					if LOCI[locus]==2:
						TOTcod2+=1
				else:
					BADcod[locus]= LOCI[locus] #used with warnings
					TOTbad+=1
					if LOCI[locus]==2:
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
	print 'N reference sequences:\t',CDSfile
	print 'N list1 (background):\t',list1
	print	'N list2 (loci selection):\t',list2
	print '- MIN required codons:\t',
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
	print 'loci with valid model:\t',len(LOCI)
	print 'of which with valid sequence:\t', Nseqs
	print 'of which passing filtering:\t',len(COUNT.keys())
	print '#loci with warnings: (NOTE: overlap is possible)'
	print '- NOT FOUND: loci without model:\t',len(member(NOTFOUND,1))
	print '- included: length not 3x:\t',len(BADlen3)
	print '- included: missing STOP:\t',len(BADstopEND)	
	print '- included: bad codons:\t',len(BADcod)
	print '- EXCLUDED: bad START:\t',len(BADstart)
	print '- EXCLUDED: internal STOP:\t',len(BADstop)
	print '- EXCLUDED: below MIN length:\t',len(BADlenMIN)
	print '#codon counts:'
	print 'N valid codons:\t',TOTcod
	print 'N bad codons - EXCLUDED:\t',TOTbad
	print
	if list2!='-':
		print '#SEQUENCE STATISTICS (LIST2)'
		print 'loci with valid model:\t',len(member(LOCI,2))
		print 'of which with valid sequence:\t',Nseqs2
		print 'of which passing filtering:\t',len(filter(lambda x:COUNT.get(x,0)!=0,member(LOCI,2)))
		print '#loci with warnings: (NOTE: overlap is possible)'
		print '- NOT FOUND: loci without model/list1 match:\t',len(member(NOTFOUND,2))
		print '- included: length not 3x:\t',len(member(BADlen3,2))
		print '- included: missing STOP:\t',len(member(BADstopEND,2))
		print '- included: bad codons:\t',len(member(BADcod,2))
		print '- EXCLUDED: bad START:\t',len(member(BADstart,2))
		print '- EXCLUDED: internal STOP:\t',len(member(BADstop,2))
		print '- EXCLUDED: below MIN length:\t',len(member(BADlenMIN,2))
		print '#codon counts:'
		print 'N valid codons:\t',TOTcod2
		print 'N bad codons - EXCLUDED:\t',TOTbad2
	
	#NORMALITY CHECKS (LIST1 ONLY)
	
	if list2=='-':
		name= '.'.join(list1.split('.')[:-1])
		if CUT>0:
			name+='.CUT'+str(CUT)
		print
		print '#NORMALITY CHECKS'
		print
		print '#effective number of codons (Nc)'
		normal(Nc,Nc.keys(),name+'.Nc')
		print
		print '#GC content'
		normal(GC,GC.keys(),name+'.GC')
		print
		print '#GC3 content'
		normal(GC3,GC3.keys(),name+'.GC3')
		if CUT==0:
			print
			print '#sequence length'
			normal(LEN,LEN.keys(),name+'.LEN')
		
	#SAMPLE MEAN VS BACKGROUND (LIST 1+2)
	
	else:
		name= '.'.join(list2.split('.')[:-1])+'_vs_'+'.'.join(list1.split('.')[:-1])
		if CUT>0:
			name+='.CUT'+str(CUT)
		print
		print '#SAMPLE vs BACKGROUND (LIST1+LIST2)'
		print
		print '#effective number of codons (Nc)'
		test(Nc,filter(lambda x:Nc.get(x,0)!=0,member(LOCI,2)),name+'.Nc','"Nc index"') 
		print
		print '#GC content'
		test(GC,filter(lambda x:GC.get(x,0)!=0,member(LOCI,2)),name+'.GC','"GC%"')
		print
		print '#GC3 content'
		test(GC3,filter(lambda x:GC3.get(x,0)!=0,member(LOCI,2)),name+'.GC3','"GC3%"')
		if CUT==0:
			print
			print '#sequence length'
			test(LEN,filter(lambda x:LEN.get(x,0)!=0,member(LOCI,2)),name+'.LEN','"sequence length (nt)"')
		print
		