#! /usr/bin/python

#V2 changes:
#- when extending upstream, stops at STOP codon (previously at START codon, causing premature interruption in case of non-start Methionine)
#- outputs new genome coordinates resulting from extraction/extension
#- adds separator (^) between upstream | within coordinates | downstream region in the output sequence (both DNA/protein)

sep= '^' #chosen because not likely to appear in sequence name, i.e. if one wants to remove it with REPLACE command in a text editor

import sys
import os
os.F_OK=1
import subprocess as sub
sub.PIPE=1

#####################################################
#TRANSLATION TABLE
#####################################################
CODE= {}
CODE['ATG']='M' #translation start
CODE['TGG']='W'
for triplet in ['TTT','TTC']:
	CODE[triplet]='F'
for triplet in ['TAT','TAC']:
	CODE[triplet]='Y'
for triplet in ['TGT','TGC']:
	CODE[triplet]='C'
for triplet in ['AAT','AAC']:
	CODE[triplet]='N'
for triplet in ['CAA','CAG']:
	CODE[triplet]='Q'
for triplet in ['GAT','GAC']:
	CODE[triplet]='D'
for triplet in ['GAA','GAG']:
	CODE[triplet]='E'
for triplet in ['AAA','AAG']:
	CODE[triplet]='K'
for triplet in ['CAT','CAC']:
	CODE[triplet]='H'
for triplet in ['ATT','ATC','ATA']:
	CODE[triplet]='I'
for triplet in ['GTT','GTC','GTA','GTG']:
	CODE[triplet]='V'
for triplet in ['ACT','ACC','ACA','ACG']:
	CODE[triplet]='T'
for triplet in ['GCT','GCC','GCA','GCG']:
	CODE[triplet]='A'
for triplet in ['CCT','CCC','CCA','CCG']:
	CODE[triplet]='P'
for triplet in ['GGT','GGC','GGA','GGG']:
	CODE[triplet]='G'
for triplet in ['TTA','TTG','CTT','CTC','CTA','CTG']:
	CODE[triplet]='L'
for triplet in ['TCT','TCC','TCA','TCG','AGT','AGC']:
	CODE[triplet]='S'
for triplet in ['CGT','CGC','CGA','CGG','AGA','AGG']:
	CODE[triplet]='R'
for triplet in ['TAA','TAG','TGA']:
	CODE[triplet]='*' #translation stop
#####################################################

compl={'A':'T','C':'G','G':'C','T':'A'} 
def revcomp(seq):
	seq2=''
	for i in range(len(seq)-1,-1,-1):
		seq2+= compl.get(seq[i].upper(),'N')
	return seq2

#####################################################

if len(sys.argv)!=5 or sys.argv[4] not in ['0','1']:
	print
	print 'usage: <sequence_file> <limits_file> <UP+seq+DW> <0/1 (EXTRACT/EXTEND)>'
	print
	print ' - UP + seq + DW extract upstream + within-limits + downstream sequence'
	print '   (UP, DW -> integer; seq -> 1/0 to include/leave out sequence within limits)'
	print ' - 0 == EXTRACT dna sequence only'
	print '   -> OUTPUT: DNA SEQUENCE'
	print ' - 1 == EXTEND protein until STOP codon found (position of 1st residue determines reading frame);'
	print '   if this is not found, EXTEND until UP/DW limits are reached'
	print '   -> OUTPUT: TRANSLATED PROTEIN SEQUENCE + CORRESPONDING DNA'
	print
	print '***IMPORTANT FOR OUTPUT FILES USAGE!***'
	print ' - remove separator ('+sep+') between upstream | within coordinates | downstream sequence'
	print '	  (both DNA/protein) before usage as regular fasta file!'
	print ' - check "CDS:NOT3x" warnings, they introduce NUMBERS (1 or 2) in the protein sequence!'
	print '   (in this case, recommended replacing translated CDS with original Blast hit sequence to account for splicing forms)'
	print
	print 'LIMITS FILE SPECIFICATIONS'
	print ' - limits should be provided as seqNAME,START,STOP,(STRAND) [string,integer,integer,(+/-)] (additional columns ignored)' 
	print '    if STRAND not given, START>STOP interpreted as reverse filament (STRAND = -)'
	print '    if given, START<STOP works with both +/-, however START>STOP will generate an error with STRAND = +'
	print ' - limits are intended as INCLUSIVE of both START:END residues, with indexing starting from 1'
	print '    (during extraction limits are converted to string notation, ex. sequence 1:4 -> string[0:4])'
	print ' - seqNAME starting with ">" tolerated'
	print ' - seqNAME with additional descriptors tolerated if following "|", ex. ATXXXX|bZIPX (bZIPX ignored during ID matching)'
	print ' - lines starting with "#" ignored'
	print
else:
	#file names
	seqfile= sys.argv[1]
	limfile= sys.argv[2]
	UgD= sys.argv[3]
	UP= abs(int(UgD.split('+')[0]))
	DW= abs(int(UgD.split('+')[2]))
	INCLUDE=False
	if UgD.split('+')[1]=='1':
		INCLUDE=True
	if sys.argv[4]=='1':
		EXTEND=True
		ntout= '.'.join(limfile.split('.')[:-1])+'_'+UgD+'_INFRAME.nt.fa'
		aaout= '.'.join(limfile.split('.')[:-1])+'_'+UgD+'_INFRAME.aa.fa'
		if os.path.isfile(ntout):
			os.remove(ntout)
		if os.path.isfile(aaout):
			os.remove(aaout)
	else:
		EXTEND=False
		outfile= '.'.join(limfile.split('.')[:-1])+'_'+UgD+'.nt.fa'
		if os.path.isfile(outfile):
			os.remove(outfile) #since outfile is written with "append", this makes sure any previous version is overwritten
	tmp= '.'.join(limfile.split('.')[:-1])+'_'+UgD+'.tmp'

	#load limits
	DB={}
	ERROR=False
	for line in open(limfile).readlines():
		x=line.split('\t')
		if len(x)>=3 and x[0][0]!='#':
			ID=x[0].split()[0]
			if ID[0]=='>':
				ID=ID[1:]
			if ID[:3]=='gb|': #to avoid problems with NCBI IDs
				ID= 'gb|'+ID.split('|')[1]
			else:
				ID= ID.split('|')[0]
			START= int(x[1])
			STOP= int(x[2])
			if len(x)>3 and x[3] in ['+','-']:
				STRAND= x[3]
				if STRAND=='+' and START>STOP:
					print '*** ERROR! check coordinates for ID',ID,'***'
					ERROR=True
					break
			else:
				if START>STOP:
					START,STOP=STOP,START
					STRAND='-'
				else:
					STRAND='+'
			coord= (START-1,STOP,STRAND)
			if coord not in DB.get(ID,[]):
				DB[ID]=DB.get(ID,[])+[coord] #correction from sequence notation -> string notation
	
	#extract relevant sequences
	if not ERROR: 
		#originally entirely done with grep, but Python couldn't handle the output from communicate() so now using intermediate temporary files.. annoying
		#Nlines= int(sub.Popen('cat '+seqfile+' | wc -l ',shell=True, stdout=1).communicate()[0])
		IDindex={}
		c,Nlines=0,0
		#N= sub.Popen('grep -n ">" '+seqfile+'',shell=True, stdout=1).communicate()[0].split('\n')[:-1] #OLD
		sub.call('grep -n ">" '+seqfile+' > '+tmp+'',shell=True)
		PREVIOUS= False
		for line in open(tmp).readlines()[:-1]: #last is an empty line
			Nlines+=1
			Nline= line.split(':')[0]		
			#if Nline.isdigit():
			if PREVIOUS:
				IDindex[ID]= (nID,int(Nline)) #start-stop of sequence
			ID= ':'.join(line.split(':')[1:]).split()[0][1:]
			if ID[:3]=='gb|':
				ID='gb|'+ID.split('|')[1]
			else:
				ID=ID.split('|')[0]
			if DB.get(ID,0)!=0:
				nID= int(Nline)
				PREVIOUS=True
			else:
				PREVIOUS=False
				if len(IDindex)==len(DB): #stop extraction if all relevant seq IDs have been localized
					break
		if PREVIOUS:
			IDindex[ID]= (nID,Nlines+1) #in case one of the IDs is the very last sequence in the fasta file
		#extract sequences with sed
		IDlist= DB.keys()
		IDlist.sort()
		BAD=[]
		IDlist2=[]
		for ID in IDlist:
			if IDindex.get(ID,0)==0:
				BAD.append(ID)
			else:
				IDlist2.append(ID)
		IDlist=IDlist2
		print
		if len(IDlist)==0:
			print '*** WARNING: no sequence ID match! ***'
			print
		else:
			for ID in IDlist:
				index,next=IDindex[ID][0],IDindex[ID][1]
				#OLD COMMAND causing Python memory to fill up with big files
				#SEQ= sub.Popen('grep -A '+str(next-nID-1)+' "'+ID+'" '+seqfile+'',shell=True, stdout=1).communicate()[0]
				#SEQ= ''.join(SEQ.split('\n')[1:-1])

				#NEW COMMAND
				sub.call('sed -n "'+str(index+1)+','+str(next-1)+' p" '+seqfile+' > '+tmp+'',shell=True)
				SEQ=''
				for line in open(tmp).readlines():
					SEQ+= line.strip()
				#extract relevant substrings for each sequence
				coord= DB[ID]
				coord.sort()
				for xy in coord:
					QUIT1=False
					QUIT2=False
					START= xy[0]
					STOP= xy[1]
					STRAND= xy[2]
					warning= ''
					seq=''
					if EXTEND:
						aa=''
					if STRAND=='+': #forward strand [+]
						if UP>0:
							if START-UP<0:
								warning+= ' UP:'+str(START)
								part= SEQ[:START]
								newSTART=0
							else:
								part= SEQ[START-UP:START]
								newSTART=START-UP
							if EXTEND:
								steps=0
								for i in range(len(part),2,-3):
									steps+=1
									triplet=part[i-3:i].upper()
									seq= triplet+seq	
									aa=CODE.get(triplet,'X')+aa
									if CODE.get(triplet,0)=='*':
										break
								newSTART=START-steps*3
							else:
								seq+=part
						if INCLUDE:
							seq+=sep
							if EXTEND:
								aa+=sep
								part= SEQ[START:STOP]
								for i in range(0,len(part)-2,3):
									triplet=part[i:i+3].upper()
									seq+= triplet
									aa+=CODE.get(triplet,'X')
								if len(part)%3!=0:
									if len(part)%3==2:
										warning+= ' CDS:NOT3x (2nt out)'
										aa+='2'
										seq+=part[-2:]
									else:
										warning+= ' CDS:NOT3x (1nt out)'
										aa+='1'
										seq+=part[-1]
							else:
								seq+= SEQ[START:STOP]
						if DW>0:
							seq+=sep
							if STOP+DW>len(SEQ):
								warning+= ' DW:'+str(len(SEQ)-STOP)
								part= SEQ[STOP:]
								newSTOP=len(SEQ)
							else:
								part= SEQ[STOP:STOP+DW]
								newSTOP=STOP+DW
							if EXTEND:
								aa+=sep
								steps=0
								for i in range(0,len(part)-2,3):
									steps+=1
									triplet=part[i:i+3].upper()
									seq+= triplet
									aa+=CODE.get(triplet,'X')
									if CODE.get(triplet,0)=='*':
										break
								newSTOP=STOP+steps*3
							else:
								seq+= part
					else: #reverse strand [-]
						if UP>0:
							if STOP+UP>len(SEQ):
								warning+= ' UP:'+str(len(SEQ)-STOP)
								part= SEQ[STOP:]
								newSTOP=len(SEQ)
							else:
								part= SEQ[STOP:STOP+UP]
								newSTOP=STOP+UP
							if EXTEND:
								part= revcomp(part)
								steps=0
								for i in range(len(part),2,-3):
									steps+=1
									triplet=part[i-3:i].upper()
									seq= triplet+seq	
									aa=CODE.get(triplet,'X')+aa
									if CODE.get(triplet,0)=='*':
										break
								newSTOP=STOP+steps*3
							else:
								seq+=part
						if INCLUDE:
							seq+=sep
							if EXTEND:
								aa+=sep
								part= revcomp(SEQ[START:STOP])
								for i in range(0,len(part)-2,3):
									triplet=part[i:i+3].upper()
									seq+= triplet
									aa+=CODE.get(triplet,'X')
								if len(part)%3!=0:
									if len(part)%3==2:
										warning+= ' CDS:NOT3x (2nt out)'
										aa+='2'
										seq+=part[-2:]
									else:	
										warning+= ' CDS:NOT3x (1nt out)'
										aa+='1'
										seq+=part[-1]					
							else:
								seq= SEQ[START:STOP]+seq
						if DW>0:
							seq+=sep
							if START-DW<0:
								warning+= ' DW:'+str(START)
								part= SEQ[:START]
								newSTART=0
							else:
								part= SEQ[START-DW:START]
								newSTART=START-DW
							if EXTEND:
								aa+=sep
								part= revcomp(part)
								steps=0
								for i in range(0,len(part)-2,3):
									steps+=1
									triplet=part[i:i+3].upper()
									seq+= triplet
									aa+=CODE.get(triplet,'X')
									if CODE.get(triplet,0)=='*':
										break
								newSTART=START-steps*3
							else:
								seq= part+seq
						if not EXTEND:
							seq= revcomp(seq)
					name= '>'+ID+' | '+str(START+1)+':'+str(STOP)+'['+STRAND+'] -> '+UgD+' = '+str(newSTART+1)+':'+str(newSTOP)
					if warning!='':
						warning= ' (WARNING:'+warning+')'
					if EXTEND:
						fnt=open(ntout,'a')
						fnt.write(name+warning+'\n')
						fnt.write(seq+'\n')
						fnt.close()		
						faa=open(aaout,'a')
						faa.write(name+warning+'\n')
						faa.write(aa+'\n')
						faa.close()	
					else:
						f=open(outfile,'a')
						f.write(name+warning+'\n')
						f.write(seq+'\n')
						f.close()		
					print name+'\t'+warning
					#print '    nt.\t'+seq
					#if EXTEND:
					#	print '    aa.\t'+aa
		if os.path.isfile(tmp):
			os.remove(tmp) #previously overwritten, needs to be removed after last iteration
		if BAD!=[]:
			print
			print '*****************************************'
			for ID in BAD:
				print 'no sequence retrieved for',ID
		print

