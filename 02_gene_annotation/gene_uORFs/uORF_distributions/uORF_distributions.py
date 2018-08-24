#!/usr/bin/python

#WHAT IS THIS PROGRAM FOR:
#given target loci and background (both as lists of loci IDs), 
#extracts distribution of several uORF features and computes p-values
#for the difference between the means or medians.
#uses data from supplementary table in von Arnim et al. 2014

#================================================

#CONTROL PANEL

modelfile=	'TAIR10_representative_gene_models.txt'
uORFfile= 	'von_Arnim_2014_TAIR10_uORFs.txt'

#================================================

#FUNCTIONS DEFINITION

import sys
import os
import subprocess as sub
sub.PIPE= 1

#----------------------------------------------------
# NORMALITY CHECK
#----------------------------------------------------

#calculates Shapiro-Wilk p-value
#generates Q-Q plot and histogram (100 bins) with superimposed normal distribution

def normal(locidict,name):

	alpha=0.05

	tmp=name+'.tmp'
	pdf1=name+'.QQplot.pdf'
	pdf2=name+'.histo.pdf'
	
	f=open(tmp,'w')
	for locus in locidict.keys():
		f.write(str(locidict[locus])+'\n')
	f.close()
	
	x=sub.Popen('Rscript test_normality_PDF.R '+tmp+' '+pdf1+' '+pdf2+' ', shell=True, stdout=1)
	out= x.communicate()[0].split() 
	os.remove(tmp)
	
	m=  float(out[0])
	sd= float(out[1])
	pval= float(out[2]) #Shapiro-Wilk test if N<5000, Kolmogorov-Smirnov otherwise

	print 'mean (sd):\t%1.2f (%1.2f)' % (m,sd)
	if len(locidict)<=5000:
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

	if list2==[]:
		return '-'
	else:

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
		
		if out==[]:
			print '#background (list1)'
			print 'N:\t\t',len(locidict.keys())
			print '#sample (list2)'
			print 'N:\t\t',len(list2)
			print '(NOT POSSIBLE TO COMPLETE TEST)'
		else:
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
	
#================================================

#START

import sys

if len(sys.argv)<2 or len(sys.argv)>3:
	print
	print 'usage: <list1> <list2 (OPTIONAL)> (both AGI loci)'
	print '- if both list1 and list2 given, 1 is the background and 2 is the target list;'
	print '  the program compares means/medians of various uORF features distributions'
	print '- if one list is given, the program checks the normality of the distributions for that list'
	print
else:
	bkg=sys.argv[1]
	infile='-'
	if len(sys.argv)==3:
		infile=sys.argv[2]
	
#================================================

#load background and target list gene models
	
	MOD={}
	for line in open(modelfile,'rU').readlines():
			x=line.split()
			if len(x)!=0 and x[0][:2]=='AT':
				locus=x[0].split('.')
				MOD[locus[0]]=locus[1]
	
	LOCI={} 		#locus:1 if only in background, :2 if also in target list (if given)
	for line in open(bkg,'rU').readlines():
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
					
	if infile!='-':
		for line in open(infile,'rU').readlines():
			x=line.split()[0]
			if len(x)>0 and x[0]!='#':
				if x[0]=='>':
					x=x[1:]
				if x[:2]=='AT':
					locus=x.split('|')[0].split('/')[0] #remove descriptors, add here other separators if needed
					if not len(locus.split('.'))>1: #gene model provided
						locus=locus+'.'+MOD[locus]
					if LOCI.get(locus,0)!=0:
						LOCI[locus]=2

	Nbkg=len(LOCI.keys())
	Nloci= len(filter(lambda x:LOCI[x]==2,LOCI))

#================================================

#extract uORF counts from data table (1 == bkg, 2 == target list)

	#loci counts
	N,Ns,Nm,Nw={},{},{},{}		#N uORFs per locus, ALL + strong/medium/weak
	Lmax,Lmean={},{}				#longest uORF length / average uORF length per locus
	OVmax,OVmean={},{}			#max uORF overlap / average uORF overlap per locus 
	L,Ls,Lm,Lw={},{},{},{}		#uORFs length, ALL + strong/medium/weak
	OV,OVs,OVm,OVw={},{},{},{}	#uORFs overlap, ALL + strong/medium/weak
	U=[]								#uORFs in target list
	
	for line in open(uORFfile,'rU').readlines()[1:]:
		x=line.split()
		locus=x[0]
		uORF=x[1]
		if LOCI.get(locus,0)!=0 and uORF!='-':
			N[locus]=N.get(locus,0.0)+1
			ovlen=int(x[6])
			strength=int(x[8])
			length=int(x[9])
			L[uORF]=length
			OV[uORF]=ovlen
			if LOCI[locus]==2:
				U.append(uORF)
			if strength==2:	#strong ATG context
				Ns[locus]=Ns.get(locus,0.0)+1
				Ls[uORF]=length
				OVs[uORF]=ovlen
			elif strength==1:	#medium ATG context
				Nm[locus]=Nm.get(locus,0.0)+1
				Lm[uORF]=length
				OVm[uORF]=ovlen
			else:					#weak ATG context
				Nw[locus]=Nw.get(locus,0.0)+1
				Lw[uORF]=length
				OVw[uORF]=ovlen
			Lmean[locus]=(Lmean.get(locus,[0,0])[0]+1.0,Lmean.get(locus,[0,0])[1]+length)
			if Lmax.get(locus,0)==0 or length>Lmax[locus]:
				Lmax[locus]=length
			OVmean[locus]=(OVmean.get(locus,[0,0])[0]+1.0,OVmean.get(locus,[0,0])[1]+ovlen)
			if OVmax.get(locus,0)==0 or ovlen>OVmax[locus]:
				OVmax[locus]=ovlen
	
	for locus in Lmean.keys():
		Lmean[locus]=Lmean[locus][1]/Lmean[locus][0]
	for locus in OVmean.keys():
		OVmean[locus]=OVmean[locus][1]/OVmean[locus][0]
	
#================================================

#RESULTS

	#sequence statistics
		
	print
	print '#INPUT FILES'
	print 'target list:\t',
	if infile !='-':
		print infile
		print 'background:\t',bkg
	else:
		print bkg
	
	#NORMALITY CHECKS (1 LIST ONLY)
	
	if infile=='-':
		name= '.'.join(bkg.split('.')[:-1])
		print
		print '#NORMALITY CHECKS'
		print
		print '*N uORFs (PER LOCUS)'
		normal(N,name+'.N')
		print '*uORFs LENGTH (LONGEST PER LOCUS)'
		normal(Lmax,name+'.Lmax')
		print '*uORFs LENGTH (AVERAGE PER LOCUS)'
		normal(Lmean,name+'.Lmean')
		print '*uORFs OVERLAP (GREATEST OVERLAP PER LOCUS)'
		normal(OVmax,name+'.OVmax')
		print '*uORFs OVERLAP (AVERAGE OVERLAP PER LOCUS)'
		normal(OVmean,name+'.OVmean')
		print '*uORFs LENGTH (ALL uORFs)'
		normal(L,name+'.L')
		print '*uORFs OVERLAP (ALL uORFs)'
		normal(OV,name+'.OV')
		print '*uORFs LENGTH (STRONG uORFs)'
		normal(Ls,name+'.Ls')
		print '*uORFs OVERLAP (STRONG uORFs)'
		normal(OVs,name+'.OVs')
		print '*uORFs LENGTH (MEDIUM uORFs)'
		normal(Lm,name+'.Lm')
		print '*uORFs OVERLAP (MEDIUM uORFs)'
		normal(OVm,name+'.OVm')
		print '*uORFs LENGTH (WEAK uORFs)'
		normal(Lw,name+'.Ls')
		print '*uORFs OVERLAP (WEAK uORFs)'
		normal(OVw,name+'.OVw')
		
	#SAMPLE MEAN/MEDIAN VS BACKGROUND (2 LISTS)
	
	else:
		name= '.'.join(infile.split('.')[:-1])+'_vs_'+'.'.join(bkg.split('.')[:-1])
		print
		print '#SAMPLE vs BACKGROUND'
		print
		#loci statistics
		print '*N uORFs (PER LOCUS)'
		test(N,filter(lambda x:N.get(x,'-')!='-',member(LOCI,2)),name+'.N','"N uORFs (all)"') 
		print
		print '*N STRONG uORFs (PER LOCUS)'
		test(Ns,filter(lambda x:Ns.get(x,'-')!='-',member(LOCI,2)),name+'.Ns','"N strong uORFs (all)"') 
		print
		print '*N MEDIUM uORFs (PER LOCUS)'
		test(Nm,filter(lambda x:Nm.get(x,'-')!='-',member(LOCI,2)),name+'.Nm','"N medium uORFs (all)"') 
		print
		print '*N WEAK uORFs (PER LOCUS)'
		test(Nw,filter(lambda x:Nw.get(x,'-')!='-',member(LOCI,2)),name+'.Nw','"N weak uORFs (all)"') 
		print				
		print '*uORFs LENGTH (LONGEST PER LOCUS)'
		test(Lmax,filter(lambda x:Lmax.get(x,'-')!='-',member(LOCI,2)),name+'.Lmax','"longest uORF length (nt)"') 
		print
		print '*uORFs LENGTH (AVERAGE PER LOCUS)'
		test(Lmean,filter(lambda x:Lmean.get(x,'-')!='-',member(LOCI,2)),name+'.Lmean','"average uORF length (nt)"') 
		print
		print '*uORFs OVERLAP (GREATEST OVERLAP PER LOCUS)'
		test(OVmax,filter(lambda x:OVmax.get(x,'-')!='-',member(LOCI,2)),name+'.OVmax','"closest uORF overlap with main ORF (-/+ nt)"') 
		print
		print '*uORFs OVERLAP (AVERAGE OVERLAP PER LOCUS)'
		test(OVmean,filter(lambda x:OVmax.get(x,'-')!='-',member(LOCI,2)),name+'.OVmean','"average uORF overlap with main ORF (-/+ nt)"') 
		print
		#uORFs statistics
		print '*uORFs LENGTH (ALL uORFs)'
		test(L,U,name+'.L','"uORF length (nt)"') 
		print
		print '*uORFs OVERLAP (ALL uORFs)'
		test(OV,U,name+'.OV','"uORF overlap with main ORF (+/- nt)"') 
		print
		print '*uORFs LENGTH (STRONG uORFs)'
		test(Ls,filter(lambda x:Ls.get(x,'-')!='-',U),name+'.Ls','"strong uORF length (nt)"') 
		print
		print '*uORFs OVERLAP (STRONG uORFs)'
		test(OVs,filter(lambda x:OVs.get(x,'-')!='-',U),name+'.OVs','"strong uORF overlap with main ORF (+/- nt)"') 
		print
		print '*uORFs LENGTH (MEDIUM uORFs)'
		test(Lm,filter(lambda x:Lm.get(x,'-')!='-',U),name+'.Lm','"medium uORF length (nt)"') 
		print
		print '*uORFs OVERLAP (MEDIUM uORFs)'
		test(OVm,filter(lambda x:OVm.get(x,'-')!='-',U),name+'.OVm','"medium uORF overlap with main ORF (+/- nt)"') 
		print
		print '*uORFs LENGTH (WEAK uORFs)'
		test(Lw,filter(lambda x:Lw.get(x,'-')!='-',U),name+'.Lw','"weak uORF length (nt)"') 
		print
		print '*uORFs OVERLAP (WEAK uORFs)'
		test(OVw,filter(lambda x:OVw.get(x,'-')!='-',U),name+'.OVw','"weak uORF overlap with main ORF (+/- nt)"') 
		print