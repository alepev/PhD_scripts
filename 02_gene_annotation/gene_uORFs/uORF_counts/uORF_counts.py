#!/usr/bin/python

#WHAT IS THIS PROGRAM FOR:
#given target loci and background (both as lists of loci IDs), 
#extracts several uORF statistics and p-values computed with: 
#- 1-tailed Fisher's exact test (if 2 categories)
#- Chi-square test for goodness-of-fit (if 3 categories)
#using data from supplementary table in von Arnim et al. 2014

#================================================

#CONTROL PANEL

modelfile=	'TAIR10_representative_gene_models.txt'
uORFfile= 	'von_Arnim_2014_TAIR10_uORFs.txt'

#================================================

#FUNCTIONS DEFINITION

import sys
import os
import subprocess as sub
sub.PIPE=1

def FoldFisher(N,K,n,k):
	#N == tot. genes in background
	#K == tot. genes in background with annotation X
	#n == genes in list
	#k == genes in list with annotation X
	if k==0 or n==0 or K==0 or N==0:
		x,pval= '-','-'
	else:
		FC=(float(k)/float(n))/(float(K)/float(N))
		if FC>=1:
			alt='greater'
		else:
			alt='less'
		y=sub.Popen('Rscript Fisher_1-tailed.R '+str(N)+' '+str(K)+' '+str(n)+' '+str(k)+' '+alt, shell=True, stdout=1)
		x= '%1.2f' % (FC)
		pval='%1.2e' % float(y.communicate()[0])
	return (x,pval)

#ex. of Fisher test for motif m:
#FISHER[m]= FoldFisher(Nbkg,Mbkg[m],Nloci,Mloci[m])

def ChiSQx3(A1,B1,C1,A2,B2,C2):
	#1 == counts per category in background (expected)
	#2 == counts per category in target list (observed) 
	if A1==0 or B1==0 or C1==0 or A2==0 or B2==0 or C2==0:
		pval= '-'
	else:
		y=sub.Popen('Rscript ChiSquare_goodness-of-fit_3x.R '+str(A1)+' '+str(B1)+' '+str(C1)+' '+str(A2)+' '+str(B2)+' '+str(C2)+' ', shell=True, stdout=1)
		pval='%1.2e' % float(y.communicate()[0])
	return pval
	
#================================================

#START

if len(sys.argv)!=3:
	print
	print 'usage: <target list> <background> (both AGI loci)'
	print
else:
	infile= sys.argv[1]
	bkg=sys.argv[2]
	out='.'.join(infile.split('.')[:-1])+'.uORFcounts.out'
	
#================================================

#load background and target list gene models
	
	MOD={}
	for line in open(modelfile,'rU').readlines():
			x=line.split()
			if len(x)!=0 and x[0][:2]=='AT':
				locus=x[0].split('.')
				MOD[locus[0]]=locus[1]
	
	LOCI={} 		#locus:1 if only in background, :2 if also in target list
	for line in open(bkg,'rU').readlines():
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
					
	for line in open(infile,'rU').readlines():
		x=line.split()[0].upper()
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
	N={}
	Ns,Nm,Nw,Nov,Nf={},{},{},{},{}
	#uORF counts
	Us1,Us2,Um1,Um2,Uw1,Uw2=0,0,0,0,0,0
	Uov1,Uov2=0,0
	Uf1,Uf2=0,0
	
	for line in open(uORFfile,'rU').readlines()[1:]:
		x=line.split()
		locus=x[0]
		uORF=x[1]
		TARGET=False
		if LOCI.get(locus,0)!=0 and uORF!='-':
			if LOCI[locus]==2:
				TARGET=True
			N[locus]=1
			frame=int(x[4])
			overlap=int(x[5])
			strength=int(x[8])
			if strength==2:	#strong ATG context
				Ns[locus]=1
				Us1+=1
				if TARGET:
					Us2+=1
			elif strength==1:	#medium ATG context
				Nm[locus]=1
				Um1+=1
				if TARGET:
					Um2+=1
			else:					#weak ATG context
				Nw[locus]=1
				Uw1+=1
				if TARGET:
					Uw2+=1
			if overlap==1:		#uORF overlaps main ORF
				Nov[locus]=1
				Uov1+=1
				if TARGET:
					Uov2+=1
			if frame==0:		#uORF is in frame with main ORF
				Nf[locus]=1
				Uf1+=1
				if TARGET:
					Uf2+=1
	
	N2= len(filter(lambda x:LOCI[x]==2,N))
	N1= len(N)
	Ns2= len(filter(lambda x:LOCI[x]==2,Ns))
	Nm2= len(filter(lambda x:LOCI[x]==2,Nm))
	Nw2= len(filter(lambda x:LOCI[x]==2,Nw))
	Ns1= len(Ns)
	Nm1= len(Nm)
	Nw1= len(Nw)
	Nov2= len(filter(lambda x:LOCI[x]==2,Nov))
	Nov1= len(Nov)
	Nf2= len(filter(lambda x:LOCI[x]==2,Nf))
	Nf1= len(Nf)
	U1= Us1+Um1+Uw1
	U2= Us2+Um2+Uw2
	
#================================================

#RESULTS

	f=open(out,'w')
	f.write('\n#INPUT FILES')
	f.write('\ntarget list:\t'+infile)
	f.write('\nbackground:\t'+bkg)
	f.write('\n\n#LOCI COUNTS')
	f.write('\n\n*LOCI WITH >= 1 uORF')
	FISHER= FoldFisher(Nbkg,N1,Nloci,N2)
	percent= '\t%4.2f\t%4.2f' % (100.0*N2/Nloci, 100.0*N1/Nbkg)
	f.write('\nTOTlist\tTOTbkg\tNlist\tNbkg\tN%list\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)')
	f.write('\n'+str(Nloci)+'\t'+str(Nbkg)+'\t'+str(N2)+'\t'+str(N1))
	f.write(percent+'\t'+str(FISHER[0])+'\t'+str(FISHER[1]))
	f.write('\n\n*uORF LOCI WITH >= 1 STRONG ATG CONTEXT uORF')
	FISHER= FoldFisher(N1,Ns1,N2,Ns2)
	percent= '\t%4.2f\t%4.2f' % (100.0*Ns2/N2, 100.0*Ns1/N1)
	f.write('\nTOTlist\tTOTbkg\tNlist\tNbkg\tN%list\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)')
	f.write('\n'+str(N2)+'\t'+str(N1)+'\t'+str(Ns2)+'\t'+str(Ns1))
	f.write(percent+'\t'+str(FISHER[0])+'\t'+str(FISHER[1]))
	f.write('\n\n*uORF LOCI WITH >= 1 MEDIUM ATG CONTEXT uORF')
	FISHER= FoldFisher(N1,Nm1,N2,Nm2)
	percent= '\t%4.2f\t%4.2f' % (100.0*Nm2/N2, 100.0*Nm1/N1)
	f.write('\nTOTlist\tTOTbkg\tNlist\tNbkg\tN%list\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)')
	f.write('\n'+str(N2)+'\t'+str(N1)+'\t'+str(Nm2)+'\t'+str(Nm1))
	f.write(percent+'\t'+str(FISHER[0])+'\t'+str(FISHER[1]))
	f.write('\n\n*uORF LOCI WITH >= 1 WEAK ATG CONTEXT uORF')
	FISHER= FoldFisher(N1,Nw1,N2,Nw2)
	percent= '\t%4.2f\t%4.2f' % (100.0*Nw2/N2, 100.0*Nw1/N1)
	f.write('\nTOTlist\tTOTbkg\tNlist\tNbkg\tN%list\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)')
	f.write('\n'+str(N2)+'\t'+str(N1)+'\t'+str(Nw2)+'\t'+str(Nw1))
	f.write(percent+'\t'+str(FISHER[0])+'\t'+str(FISHER[1]))
	f.write('\n\n*uORF LOCI WITH OVERLAP TO MAIN ORF')
	FISHER= FoldFisher(N1,Nov1,N2,Nov2)
	percent= '\t%4.2f\t%4.2f' % (100.0*Nov2/N2, 100.0*Nov1/N1)
	f.write('\nTOTlist\tTOTbkg\tNlist\tNbkg\tN%list\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)')
	f.write('\n'+str(N2)+'\t'+str(N1)+'\t'+str(Nov2)+'\t'+str(Nov1))
	f.write(percent+'\t'+str(FISHER[0])+'\t'+str(FISHER[1]))
	f.write('\n\n*uORF LOCI IN FRAME WITH MAIN ORF')
	FISHER= FoldFisher(N1,Nf1,N2,Nf2)
	percent= '\t%4.2f\t%4.2f' % (100.0*Nf2/N2, 100.0*Nf1/N1)
	f.write('\nTOTlist\tTOTbkg\tNlist\tNbkg\tN%list\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)')
	f.write('\n'+str(N2)+'\t'+str(N1)+'\t'+str(Nf2)+'\t'+str(Nf1))
	f.write(percent+'\t'+str(FISHER[0])+'\t'+str(FISHER[1]))
	f.write('\n\n*uORF LOCI ATG CONTEXT STRENGTH (STRONGEST uORF)')
	f.write('\nLIST:strong\tmedium\tweak\tBKG:STRONG\tMEDIUM\tWEAK\tp-value(Chi-square)')
	f.write('\n'+str(Ns2)+'\t'+str(Nm2)+'\t'+str(Nw2)+'\t'+str(Ns1)+'\t'+str(Nm1)+'\t'+str(Nw1))
	f.write('\t'+str(ChiSQx3(Ns1,Nm1,Nw1,Ns2,Nm2,Nw2)))
	
	f.write('\n\n\n#uORF COUNTS')
	f.write('\n\n*uORFs WITH STRONG ATG CONTEXT')
	FISHER= FoldFisher(U1,Us1,U2,Us2)
	percent= '\t%4.2f\t%4.2f' % (100.0*Us2/U2, 100.0*Us1/U1)
	f.write('\nTOTlist\tTOTbkg\tNlist\tNbkg\tN%list\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)')
	f.write('\n'+str(U2)+'\t'+str(U1)+'\t'+str(Us2)+'\t'+str(Us1))
	f.write(percent+'\t'+str(FISHER[0])+'\t'+str(FISHER[1]))
	f.write('\n\n*uORFs WITH MEDIUM ATG CONTEXT')
	FISHER= FoldFisher(U1,Um1,U2,Um2)
	percent= '\t%4.2f\t%4.2f' % (100.0*Um2/U2, 100.0*Um1/U1)
	f.write('\nTOTlist\tTOTbkg\tNlist\tNbkg\tN%list\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)')
	f.write('\n'+str(U2)+'\t'+str(U1)+'\t'+str(Um2)+'\t'+str(Um1))
	f.write(percent+'\t'+str(FISHER[0])+'\t'+str(FISHER[1]))
	f.write('\n\n*uORFs WITH WEAK ATG CONTEXT')
	FISHER= FoldFisher(U1,Uw1,U2,Uw2)
	percent= '\t%4.2f\t%4.2f' % (100.0*Uw2/U2, 100.0*Uw1/U1)
	f.write('\nTOTlist\tTOTbkg\tNlist\tNbkg\tN%list\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)')
	f.write('\n'+str(U2)+'\t'+str(U1)+'\t'+str(Uw2)+'\t'+str(Uw1))
	f.write(percent+'\t'+str(FISHER[0])+'\t'+str(FISHER[1]))
	f.write('\n\n*uORFs OVERLAPPING MAIN ORF')
	FISHER= FoldFisher(U1,Uov1,U2,Uov2)
	percent= '\t%4.2f\t%4.2f' % (100.0*Uov2/U2, 100.0*Uov1/U1)
	f.write('\nTOTlist\tTOTbkg\tNlist\tNbkg\tN%list\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)')
	f.write('\n'+str(U2)+'\t'+str(U1)+'\t'+str(Uov2)+'\t'+str(Uov1))
	f.write(percent+'\t'+str(FISHER[0])+'\t'+str(FISHER[1]))
	f.write('\n\n*uORFs IN FRAME WITH MAIN ORF')
	FISHER= FoldFisher(U1,Us1,U2,Us2)
	percent= '\t%4.2f\t%4.2f' % (100.0*Uf2/U2, 100.0*Uf1/U1)
	f.write('\nTOTlist\tTOTbkg\tNlist\tNbkg\tN%list\tN%bkg\tfoldEnrich\tp-value(Fisher_1-tailed)')
	f.write('\n'+str(U2)+'\t'+str(U1)+'\t'+str(Uf2)+'\t'+str(Uf1))
	f.write(percent+'\t'+str(FISHER[0])+'\t'+str(FISHER[1]))
	f.write('\n\n*uORFs ATG CONTEXT STRENGTH')
	f.write('\nLIST:strong\tmedium\tweak\tBKG:STRONG\tMEDIUM\tWEAK\tp-value(Chi-square)')
	f.write('\n'+str(Us2)+'\t'+str(Um2)+'\t'+str(Uw2)+'\t'+str(Us1)+'\t'+str(Um1)+'\t'+str(Uw1))
	f.write('\t'+str(ChiSQx3(Us1,Um1,Uw1,Us2,Um2,Uw2)))
	f.close()
