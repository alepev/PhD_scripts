#! /usr/bin/python

import sys
import subprocess as sub
import os
os.F_OK=1

#RUNS MEME ON ALL THE INPUT FILES CONTAINED IN A THE SPECIFIED INPUT FOLDER

#IMPORTANT: USEFUL WHEN MULTIPLE SHORT RUNS OF MEME, OTHERWISE BETTER INDIVIDUAL JOBS STARTED ON THE SAME MUTANT (1 PROCESS PER CORE)
#(parallel version currently not supported due to missing job management system on mutants)

MINseq=10 #IMPORTANT: INDICATES MIN. SEQUENCES FOR LIST TO BE PROCESSED!

MEMEdir= '/home/alessia/02_PROJECT_ANALYSES/06_MOTIF_SEARCH/01_MEME/' #bfiles should be contained in here!
outdir= 'RESULTS_UNIQUEeOVERLAP/' #change name according to your needs; folder will fall under MEMEdir

if len(sys.argv)!=3:
	print 'usage: <input folder> <background file>'
else:
	indir=sys.argv[1]
	bck=sys.argv[2] #background file, ex. background_500up_order2.txt
	inlist=os.listdir(indir)

	opslist=[]
	#MEME INPUT PARAMETERS - see below to match position in the list with actual function
	opslist.append(['zoops',10,'',10,1000,6,10,True,True,False,bck])
	opslist.append(['anr',10,'',10,1000,6,10,True,True,False,bck])

	for op in opslist:
	#PARAMETERS ASSIGNMENT
		mod=op[0] 		#[oops|zoops|anr]
		nmotifs=str(op[1])
		evt=str(op[2])		#e-value threshold (default infinite)
		minsites=str(op[3])
		maxsites=str(op[4])
		minw=str(op[5])
		maxw=str(op[6])
		dna=op[7]
		revcomp=op[8]
		pal=op[9]
		bfile=op[10]
		######################
	
		if indir[-1]=='/':
			indir=indir[:-1]
		outsubdir=MEMEdir+outdir+indir.split('/')[-1]+'_'+mod+'_Nsites'+minsites
		if maxsites!='':
			outsubdir+='-'+maxsites
			maxsites=' -maxsites '+maxsites
		outsubdir+='_width'+minw+'-'+maxw

		if dna:
			dna=' -dna'
			if revcomp:
				revcomp= ' -revcomp'
				outsubdir+='_bothStrands'
			else:	
				revcomp=''	
			if pal:
				pal=' -pal'
				outsubdir+='_palindrome'
			else:
				pal=''
		else:
			dna=''
			revcomp=''
			pal=''
	
		if evt!='':
			evt=' -evt '+evt
		if bfile!='':
			bfile=' -bfile '+bfile
			outsubdir+='_bkg-corrected'
		
		outsubdir+='/'
		if not os.access(outsubdir,1):
			os.mkdir(outsubdir)	
	
		for dataset in inlist:
			c=0
			for line in open(indir+'/'+dataset).readlines():
				x=line.split()
				if len(x)>0 and x[0][0]=='>':
					c+=1
			if c < MINseq:
				print '\n***file',dataset,'skipped, too few sequences found (MIN.',MINseq,'required)***\n'
				continue
			filefolder=outsubdir+dataset[:-4]
			args='meme '+indir+'/'+dataset+dna+revcomp+pal+' -maxsize 1000000 -mod '+mod+' -nmotifs '+nmotifs+' -minsites '+minsites+maxsites+' -minw '+minw+' -maxw '+maxw+bfile+evt
			print args
			args+=' -oc '+filefolder
			sub.call(args,shell=True)
			print '...done'
