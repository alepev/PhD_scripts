#!/usr/bin/python

#NOTE - modified to:
#- convert loci to uppercase
#- GC content as %, not ratio (i.e. 2 additional informative values)

import sys
import os
import subprocess as sub
sub.PIPE= 1

MIN=0		#generate counts only for genes with N nts >= MIN

if len(sys.argv) <4 or len(sys.argv)>5:
	print
	print 'usage: <DNA seqfile (ex. UTRs)> <list 1> <list 2> <background (OPTIONAL)>'
	print 'compares LENGTH and GC% of the two lists as independent samples'
	print '(NOTE: sequences do not need to be ORFs/CDSs!)'
	print '(both parametric (Z-test,T-test) and non-parametric (Wilcoxon) tests are used;'
	print 'if provided, background is used to compute population variance for the Z-test)'
	print
	if MIN>0:
		print
		print 'NOTE:\nsequences with less than',MIN,'nts are excluded from the analysis (see code to modify this behavior)'
	print
	
else:

	modelfile=	'TAIR10_representative_gene_models.txt'

	seqfile= sys.argv[1]
	list1= sys.argv[2] 
	list2= sys.argv[3]
	bkg='-'
	if len(sys.argv)==5:
		bkg=sys.argv[4]
	
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
 		print out
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
		
	#generate list of keys for which dict[key] <= or >= value
	def membercfr(dict,value,cfr):
		if cfr=='<=':
			return filter(lambda x:dict[x]<=value,dict)
		elif cfr=='>=':
			return filter(lambda x:dict[x]>=value,dict)
				
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

	#scan sequence file for length & GC counts

	GC={}
	LEN={}	
	BADnt={}			#included with warnings
	BADlenMIN={}	#excluded
	NTtot1,NTtot2=0,0
	BADtot1,BADtot2=0,0
	Nseqs1,Nseqs2=0,0
	locus=''
	seq=''
	READ=False
	for line in open(seqfile,'rU').readlines():
			if line[0]=='>':
				if locus!='' and READ: #if new ID, check previous sequence
					if seq[-1]=='*':
						seq=seq[:-1]
					if len(seq)<MIN and LOCI.get(locus,-1)!=-1:
						BADlenMIN[locus]=LOCI[locus]
					else:
						l=len(seq)
						gc,bad=0,0
						for i in range(0,l):
							if seq[i] in ['G','C']:
								gc+=1
							elif seq[i] not in ['A','T']: #bad nt
								bad+=1	
								if LOCI.get(locus,-1)!=-1:
									BADnt[locus]=LOCI[locus]
						if LOCI.get(locus,-1)!=-1:
							if LOCI[locus]<=1:
								NTtot1+=l
								BADtot1+=bad
							if LOCI[locus]>=1:
								NTtot2+=l
								BADtot2+=bad
						LEN[locus]=l
						GC[locus]=gc*100.0/(l-bad)		
										
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
		if seq[-1]=='*':
			seq=seq[:-1]
		if len(seq)<MIN and LOCI.get(locus,-1)!=-1:
			BADlenMIN[locus]=LOCI[locus]
		else:
			l=len(seq)
			gc,bad=0.0,0.0
			for i in range(0,l):
				if seq[i] in ['G','C']:
					gc+=1
				elif seq[i] not in ['A','T']: #bad nt
					bad+=1	
					if LOCI.get(locus,-1)!=-1:
						BADnt[locus]=LOCI[locus]
			if LOCI.get(locus,-1)!=-1:
				if LOCI[locus]<=1:
					NTtot1+=l
					BADtot1+=bad
				if LOCI[locus]>=1:
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
	print 'reference sequences:\t',seqfile
	print 'list1:\t',list1
	print	'list2:\t',list2
	print '- background population (Z-test variance):\t',bkg
	print '- MIN required nts:\t',
	if MIN>0:
		print MIN
	else:
		print '-'
	print
	print '#SEQUENCE STATISTICS (LIST1)'
	print 'loci with valid model:\t',len(membercfr(LOCI,1,'<='))
	print 'of which with valid sequence:\t', Nseqs1
	print 'of which passing filtering:\t',len(filter(lambda x:LEN.get(x,0)!=0,membercfr(LOCI,1,'<=')))
	print '#loci with warnings:'
	print '- included: contain invalid nts:\t',len(membercfr(BADnt,1,'<='))
	print '- EXCLUDED: below MIN length:\t',len(membercfr(BADlenMIN,1,'<='))
	print '#nt counts:'
	print 'N valid nts:\t',NTtot1
	print 'N bad nts (EXCLUDED):\t',BADtot1
	print
	print '#SEQUENCE STATISTICS (LIST2)'
	print 'loci with valid model:\t',len(membercfr(LOCI,1,'>='))
	print 'of which with valid sequence:\t', Nseqs2
	print 'of which passing filtering:\t',len(filter(lambda x:LEN.get(x,0)!=0,membercfr(LOCI,1,'>=')))
	print '#loci with warnings:'
	print '- included: contain invalid nts:\t',len(membercfr(BADnt,1,'>='))
	print '- EXCLUDED: below MIN length:\t',len(membercfr(BADlenMIN,1,'>='))
	print '#nt counts:'
	print 'N valid nts:\t',NTtot2
	print 'N bad nts (EXCLUDED):\t',BADtot2
	print
		
	#LIST1 VS LIST2

	name= '.'.join(list1.split('.')[:-1])+'_vs_'+'.'.join(list2.split('.')[:-1])
	print
	print '#LIST1 VS LIST2'
	print
	print '#SEQUENCE LENGTH'
	print
	test(LEN,LOCI,name+'.LEN','"sequence length (nt)"',BKG) 
	print
	print '#GC CONTENT'
	print
	test(GC,LOCI,name+'.GC','"GC%"',BKG)
	print