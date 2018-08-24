#!/usr/bin/python

#NOTE - modified to:
#- convert loci to uppercase
#- NOT remove temporary lists when single file (for re-openining and plotting in R)
#- GC content as %, not ratio (i.e. 2 additional informative values)

import sys
import os
import subprocess as sub
sub.PIPE= 1

MIN=0		#generate counts only for genes with N nts >= MIN

if len(sys.argv)<3 or len(sys.argv)>4:
	print
	print 'usage: <DNA seqfile (ex. UTRs)> <list 1> <list 2 (OPTIONAL)>'
	print 'extracts length and GC content for each gene (sequences do not need to be ORFs!);'
	print '- if only list1, tests for normality of the distributions (qqplot, histogram, and Shapiro normality test)'
	print '- if both list1 and list2, assumes list1 as population and tests for significantly different mean/median of list2'
	print '(both parametric (Z-test,T-test) and non-parametric (Wilcoxon) tests are used,'
	print ' see output NORMALITY CHECK plots to decide what test result to look at for each gene feature)'
	print
	if MIN>0:
		print
		print 'NOTE:\nsequences with less than',MIN,'nts are excluded from the analysis (see code to modify this behavior)'
	print
	
else:

	modelfile=	'TAIR10_representative_gene_models.txt'

	seqfile= sys.argv[1]
	list1= sys.argv[2] #background, ex. genes represented on the microarray
	list2='-'
	if len(sys.argv)==4:
		list2=sys.argv[3]
	
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
#		os.remove(tmp)
		
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
		print
	
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
		print
		
	#----------------------------------------------------
	# GENERAL FUNCTIONS
	#----------------------------------------------------
	
	#generate list of keys for which dict[key]==value
	def member(dict,value):
		return filter(lambda x:dict[x]==value,dict)
				
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
		
	#scan sequence file for length & GC counts

	GC={}
	LEN={}	
	BADnt={}			#included with warnings
	BADlenMIN={}	#excluded
	NTtot,NTtot2=0,0
	BADtot,BADtot2=0,0
	Nseqs,Nseqs2=0,0
	locus=''
	seq=''
	READ=False
	for line in open(seqfile,'rU').readlines():
			if line[0]=='>':
				if locus!='' and READ: #if new ID, check previous sequence
					if seq[-1]=='*':
						seq=seq[:-1]
					if len(seq)<MIN:
						BADlenMIN[locus]=LOCI[locus]
					else:
						l=len(seq)
						gc,bad=0,0
						for i in range(0,l):
							if seq[i] in ['G','C']:
								gc+=1
							elif seq[i] not in ['A','T']: #bad nt
								bad+=1	
								BADnt[locus]=LOCI[locus]
						NTtot+=l
						BADtot+=bad
						if LOCI[locus]==2:
							NTtot2+=l
							BADtot2+=bad
						LEN[locus]=l
						GC[locus]=gc*100.0/(l-bad)						
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
		if seq[-1]=='*':
			seq=seq[:-1]
		if len(seq)<MIN:
			BADlenMIN[locus]=LOCI[locus]
		else:
			l=len(seq)
			gc,bad=0,0
			for i in range(0,l):
				if seq[i] in ['G','C']:
					gc+=1
				elif seq[i] not in ['A','T']: #bad nt
					bad+=1	
					BADnt[locus]=LOCI[locus]
			NTtot+=l
			BADtot+=bad
			if LOCI[locus]==2:
				NTtot2+=l
				BADtot2+=bad
			LEN[locus]=l
			GC[locus]=gc*100.0/(l-bad)						
							
	#####################################################
	# OUTPUT
	#####################################################
				
	#sequence statistics
		
	print
	print '#INPUT PARAMETERS'
	print 'N reference sequences:\t',seqfile
	print 'N list1 (background):\t',list1
	print	'N list2 (loci selection):\t',list2
	print '- MIN required nts:\t',
	if MIN>0:
		print MIN
	else:
		print '-'
	print
	print '#SEQUENCE STATISTICS (LIST1)'
	print 'loci with valid model:\t',len(LOCI)
	print 'of which with valid sequence:\t', Nseqs
	print 'of which passing filtering:\t',len(LEN.keys())
	print '#loci with warnings:'
	print '- included: contain invalid nts:\t',len(BADnt)
	print '- EXCLUDED: below MIN length:\t',len(BADlenMIN)
	print '#nt counts:'
	print 'N valid nts:\t',NTtot
	print 'N bad nts (EXCLUDED):\t',BADtot
	print
	if list2!='-':
		print '#SEQUENCE STATISTICS (LIST2)'
		print 'loci with valid model:\t',len(member(LOCI,2))
		print 'of which with valid sequence:\t', Nseqs2
		print 'of which passing filtering:\t',len(filter(lambda x:LEN.get(x,0)!=0,member(LOCI,2)))
		print '#loci with warnings:'
		print '- included: contain invalid nts:\t',len(member(BADnt,2))
		print '- EXCLUDED: below MIN length:\t',len(member(BADlenMIN,2))
		print '#nt counts:'
		print 'N valid nts:\t',NTtot2
		print 'N bad nts (EXCLUDED):\t',BADtot2
		print
	
	#NORMALITY CHECKS (LIST1 ONLY)
	
	if list2=='-':
		name= '.'.join(list1.split('.')[:-1])
		print
		print '#NORMALITY CHECKS'
		print
		print '#SEQUENCE LENGTH'
		normal(LEN,LEN.keys(),name+'.LEN')
		print '#GC CONTENT'
		normal(GC,GC.keys(),name+'.GC')
		print
		
	#SAMPLE MEAN VS BACKGROUND (LIST 1+2)
	
	else:
		name= '.'.join(list2.split('.')[:-1])+'_vs_'+'.'.join(list1.split('.')[:-1])
		print
		print '#SAMPLE vs BACKGROUND (LIST1+LIST2)'
		print
		print '#SEQUENCE LENGTH'
		print
		test(LEN,filter(lambda x:LEN.get(x,0)!=0,member(LOCI,2)),name+'.LEN','"sequence length (nt)"') 
		print
		print '#GC CONTENT'
		print
		test(GC,filter(lambda x:GC.get(x,0)!=0,member(LOCI,2)),name+'.GC','"GC%"')
		print
		