#! /usr/bin/python

#v2.0: fixed a problem with sequences starting and/or ending with '*' with option INCLUDE on, which would cause their exclusion from the output

#NOTE: after bug fix there was a change in warnings (plus retrieved shorter seqs), perhaps is good but check just in case

import sys
import os
os.F_OK=1
import subprocess as sub
sub.PIPE=1

#####################################################
#TRANSLATION TABLE (STANDARD CODE)
#####################################################
CODE= {}
CODE['ATG']='M' #translation start
CODE['TGG']='W'
for coding in ['TTT','TTC']:
	CODE[coding]='F'
for coding in ['TAT','TAC']:
	CODE[coding]='Y'
for coding in ['TGT','TGC']:
	CODE[coding]='C'
for coding in ['AAT','AAC']:
	CODE[coding]='N'
for coding in ['CAA','CAG']:
	CODE[coding]='Q'
for coding in ['GAT','GAC']:
	CODE[coding]='D'
for coding in ['GAA','GAG']:
	CODE[coding]='E'
for coding in ['AAA','AAG']:
	CODE[coding]='K'
for coding in ['CAT','CAC']:
	CODE[coding]='H'
for coding in ['ATT','ATC','ATA']:
	CODE[coding]='I'
for coding in ['GTT','GTC','GTA','GTG']:
	CODE[coding]='V'
for coding in ['ACT','ACC','ACA','ACG']:
	CODE[coding]='T'
for coding in ['GCT','GCC','GCA','GCG']:
	CODE[coding]='A'
for coding in ['CCT','CCC','CCA','CCG']:
	CODE[coding]='P'
for coding in ['GGT','GGC','GGA','GGG']:
	CODE[coding]='G'
for coding in ['TTA','TTG','CTT','CTC','CTA','CTG']:
	CODE[coding]='L'
for coding in ['TCT','TCC','TCA','TCG','AGT','AGC']:
	CODE[coding]='S'
for coding in ['CGT','CGC','CGA','CGG','AGA','AGG']:
	CODE[coding]='R'
for coding in ['TAA','TAG','TGA']:
	CODE[coding]='*' #translation stop
#####################################################

compl={'A':'T','C':'G','G':'C','T':'A'} 
def revcomp(seq):
	seq2=''
	for i in range(len(seq)-1,-1,-1):
		seq2+= compl.get(seq[i].upper(),'N')
	return seq2

#####################################################

if len(sys.argv)not in [5,6] or sys.argv[4] not in ['0','1']:
	print
	print 'usage: <DNA sequence_file> <limits_file> <UP+seq+DW> <0/1 (EXTRACT/EXTEND)> <PROTEIN subsequence_file (optional)>'
	print
	print ' * UP + seq + DW extract upstream + within-limits + downstream sequence'
	print '   (UP, DW -> integer; seq -> 1/0 to include/leave out sequence within limits)'
	print
	print '   optional PROTEIN subsequence file can be provided for use with seq==1 in UP+seq+DW (UP and/or DW != 0) and EXTEND == 1 - with other settings it will be ignored'
	print '   (protein sequence from file is used in place of translated within-limits dna sequence, and joined with up/downstream translated sequence)'
	print
	print ' * mode:'
	print ' - 0 == EXTRACT DNA sequence according to provided limits'
	print '   -> OUTPUT: DNA SEQUENCE'
	print ' - 1 == EXTEND protein in frame until Stop codon or UP,DW limits are reached'
	print '   (START in limits file determines reading frame)'
	print '   -> OUTPUT: TRANSLATED PROTEIN SEQUENCE + CORRESPONDING DNA'
	print '   * if PROTEIN subs. file provided: EXTEND protein in frame from it'
	print
	print '***LIMITS FILE SPECIFICATIONS***'
	print ' - limits should be provided as seqNAME,START,STOP,(STRAND) [string,integer,integer,(+/-)] (additional columns ignored)' 
	print '    if STRAND not given, START>STOP interpreted as reverse filament (STRAND = -)'
	print '    if given, START<STOP works with both +/-, however START>STOP will generate an error with STRAND = +'
	print ' - limits are intended as INCLUSIVE of both START:END residues, with indexing starting from 1'
	print '    (during extraction limits are converted to string notation, ex. sequence 1:4 -> string[0:4])'
	print ' - seqNAME starting with ">" tolerated'
	print ' - seqNAME with additional descriptors tolerated if following "|", ex. ATXXXX|bZIPX (bZIPX ignored during ID matching)'
	print ' - lines starting with "#" ignored'
	print
	print '***(OPTIONAL) PROTEIN FILE SPECIFICATIONS***'
#	print ' - file extension should be ".aa.fa" to be properly recognized as protein sequence file'
	print ' - protein subsequence IDs and limits should match those provided in limits file' 
	print
else:
	#file names
	seqfile= sys.argv[1]
	limfile= sys.argv[2]
	UgD= sys.argv[3]
	UP= int(UgD.split('+')[0])
	DW= int(UgD.split('+')[2])
	INCLUDE=False
	if UgD.split('+')[1]=='1':
		INCLUDE=True
	if sys.argv[4]=='1':
		EXTEND=True
		USEFILE=False
		if len(sys.argv)==6:
			protfile= sys.argv[5]
			USEFILE=True
		if USEFILE:
			ntout= '.'.join(limfile.split('.')[:-1])+'_'+UgD+'_EXTENDED_FILE.nt.fa'
			aaout= '.'.join(limfile.split('.')[:-1])+'_'+UgD+'_EXTENDED_FILE.aa.fa'
		else:
			ntout= '.'.join(limfile.split('.')[:-1])+'_'+UgD+'_EXTENDED.nt.fa'
			aaout= '.'.join(limfile.split('.')[:-1])+'_'+UgD+'_EXTENDED.aa.fa'
		if os.path.isfile(ntout):
			os.remove(ntout)
		if os.path.isfile(aaout):
			os.remove(aaout)
	else:
		EXTEND=False
		outfile= '.'.join(limfile.split('.')[:-1])+'_'+UgD+'.nt.fa'
		if os.path.isfile(outfile):
			os.remove(outfile) #since outfile is written with "append", this makes sure any previous version is overwritten

	if USEFILE and (not INCLUDE or (UP==0 and DW==0)):# or protfile[-6:]!='.aa.fa'):
		print
		print '*** ERROR: check settings/optional file extension (not possible to use '+protfile.split('/')[-1]+' for analysis) ***'
		print
	else:
		tmp= '.'.join(limfile.split('.')[:-1])+'_'+UgD+'.tmp'
		N=0
		#load limits
		DB={}
		ERROR=False
		for line in open(limfile).readlines():
			x=line.split('\t')
			if len(x)>=3 and x[0][0]!='#':
				ID=x[0].split()[0]
				if ID[0]=='>':
					ID=ID[1:]
				#if ID[:3]=='gb|': #to avoid problems with NCBI IDs
				#	ID= 'gb|'+ID.split('|')[1]
				#else:
				#	ID= ID.split('|')[0]
				START= int(x[1])
				STOP= int(x[2])
				if len(x)>3 and x[3] in ['+','-']:
					STRAND= x[3]
					if STRAND=='+' and START>STOP:
						print '*** ERROR! check coordinates for ID',ID,'***'
						ERROR=True
						break
				else: #if no strand info, inferred from START,STOP limits
					if START>STOP:
						START,STOP=STOP,START
						STRAND='-'
					else:
						STRAND='+'
				coord= (START-1,STOP,STRAND) #correction from sequence notation -> string notation
				if coord not in DB.get(ID,[]):
					DB[ID]=DB.get(ID,[])+[coord] 
					N+=1
					
		#extract relevant sequences
		if not ERROR: 
			if USEFILE:
				FILEseq={}
				seq=''
				READ= False
				for line in open(protfile).readlines():
					if len(line.split())>0:
						if line[0]=='>':
							if READ:
								FILEseq[(ID,coord)]=seq
							seq=''
							READ=False
							ID= line.split()[0][1:] #' | '.join(line[1:-1].split(' | ')[:-1]) #removes ">" and "\n", extracts sequence name which matches limits file
							if DB.get(ID,0)!=0:
								coord= line.split(' | ')[-1] #extracts coordinates which match limits file
								start=int(coord.split(':')[0])
								stop= int(coord.split(':')[1].split()[0])
								strand=coord.split()[1][1] #removes square brackets
								coord= (start-1,stop,strand) #correction from sequence notation -> string notation
								if coord in DB[ID] and FILEseq.get((ID,coord),0)==0:
									READ=True
						elif READ:
							seq+=line.split()[0]
				if READ:
					FILEseq[(ID,coord)]=seq	
					
			#DISMISSED: originally done with grep, but Python couldn't handle the output from communicate() - now using intermediate temporary files
			#(DISMISSED) Nlines= int(sub.Popen('cat '+seqfile+' | wc -l ',shell=True, stdout=1).communicate()[0])
			IDindex={}
			c,Nlines=0,0
			#(DISMISSED) N= sub.Popen('grep -n ">" '+seqfile+'',shell=True, stdout=1).communicate()[0].split('\n')[:-1]
			sub.call('grep -n ">" '+seqfile+' > '+tmp+'',shell=True)
			PREVIOUS= False
			for line in open(tmp).readlines():
				Nlines+=1
				Nline= line.split(':')[0]		
				if PREVIOUS:
					IDindex[ID]= (nID,int(Nline)) #start-stop of sequence
				ID= ':'.join(line.split(':')[1:]).split()[0][1:]#[1:-1]
				#if ID[:3]=='gb|': #to avoid problems with NCBI IDs
				#	ID='gb|'+ID.split('|')[1]
				#else:
				#	ID=ID.split('|')[0]
				if DB.get(ID,0)!=0:
					nID= int(Nline)
					PREVIOUS=True
				else:
					PREVIOUS=False
					if len(IDindex)==len(DB): #stop extraction if all relevant seq IDs have been localized
						break
			if PREVIOUS:
				IDindex[ID]= (nID,Nlines+1) #in case one of the IDs is the very last sequence in the fasta file
			#check if sequences in limits file actually present in seqfile
			IDlist= DB.keys()
			IDlist.sort()
			BAD,BAD2=[],[]
			IDlist2=[]
			for ID in IDlist:
				if IDindex.get(ID,0)==0:
					BAD.append(ID)
				else:
					IDlist2.append(ID)
			IDlist=IDlist2
			#check if sequences in limits file actually present in provided protfile
			if USEFILE:
				DB2={}
				for ID in IDlist:
					for coord in DB[ID]:
						if FILEseq.get((ID,coord),0)==0:
							BAD2.append((ID,coord))
						else:
							DB2[ID]= DB2.get(ID,[])+[coord]
				DB=DB2
				IDlist=DB.keys()
			print
			#extract sequences with sed
			if len(IDlist)==0:
				print '*** WARNING: no sequence ID match! ***'
				print
			else:
				for ID in IDlist:
					#print ID
					index,next=IDindex[ID][0],IDindex[ID][1]

					#DISMISSED: command causing Python memory to fill up with big files
					#(DISMISSED) SEQ= sub.Popen('grep -A '+str(next-nID-1)+' "'+ID+'" '+seqfile+'',shell=True, stdout=1).communicate()[0]
					#(DISMISSED) SEQ= ''.join(SEQ.split('\n')[1:-1])

					#NEW COMMAND:
					sub.call('sed -n "'+str(index+1)+','+str(next-1)+' p" '+seqfile+' > '+tmp+'',shell=True)
					SEQ=''
					for line in open(tmp).readlines():
						SEQ+= line.strip()

					#extract relevant substrings for each sequence
					coord= DB[ID]
					coord.sort()
					for xy in coord:
						START= xy[0]
						STOP= xy[1]
						STRAND= xy[2]
						warning= ''
						seq=''
						if EXTEND:
							aa=''
						if STRAND=='+': #forward strand [+]
							if UP>0 and not ((USEFILE and FILEseq[(ID,(START,STOP,STRAND))][0]=='*') or (not USEFILE and CODE.get(SEQ[START:START+3],0)=='*')):
								if START-UP<0:
									warning+= ' UP:'+str(START)
									part= SEQ[:START]
								else:
									part= SEQ[START-UP:START]
								if EXTEND:								
									s1= START+1
									for i in range(len(part),2,-3):
										s1-=3
										coding=part[i-3:i].upper()
										seq= coding+seq	
										aa=CODE.get(coding,'X')+aa
										if CODE.get(coding,0)=='*':
											break
								else:
									seq+=part
							else:
								s1=START+1
							if INCLUDE:
								if USEFILE:
									seq+= SEQ[START:STOP]
									aa+= FILEseq[(ID,(START,STOP,STRAND))]
								elif EXTEND:
									part= SEQ[START:STOP]
									for i in range(0,len(part)-2,3):
										coding=part[i:i+3].upper()
										seq+= coding
										aa+=CODE.get(coding,'X')
									if len(part)%3!=0:
										warning+= ' CDS:NOT3x'
										if len(part)%3==2:
											aa+='2'
											seq+=part[-2:]
										else:
											aa+='1'
											seq+=part[-1]
								else:
									seq+= SEQ[START:STOP]
							if DW>0 and not ((USEFILE and FILEseq[(ID,(START,STOP,STRAND))][-1]=='*') or (not USEFILE and CODE.get(SEQ[STOP-3:STOP],0)=='*')):
								if STOP+DW>len(SEQ):
									warning+= ' DW:'+str(len(SEQ)-STOP)
									part= SEQ[STOP:]
								else:
									part= SEQ[STOP:STOP+DW]
								if EXTEND:
									s2=STOP
									for i in range(0,len(part)-2,3):
										s2+=3
										coding=part[i:i+3].upper()
										seq+= coding
										aa+=CODE.get(coding,'X')
										if CODE.get(coding,0)=='*':
											break
								else:
									seq+= part
							else:
								s2=STOP
						else: #reverse strand [-]
							if UP>0 and not ((USEFILE and FILEseq[(ID,(START,STOP,STRAND))][0]=='*') or (not USEFILE and CODE.get(revcomp(SEQ[STOP-3:STOP]),0)=='*')):
								if STOP+UP>len(SEQ):
									warning+= ' UP:'+str(len(SEQ)-STOP)
									part= SEQ[STOP:]
								else:
									part= SEQ[STOP:STOP+UP]
								if EXTEND:
									part= revcomp(part)
									s2=STOP
									for i in range(len(part),2,-3):
										s2+=3
										coding=part[i-3:i].upper()
										seq= coding+seq	
										aa=CODE.get(coding,'X')+aa
										if CODE.get(coding,0)=='*':
											break
								else:
									seq+=part
							else:
								s2=STOP
							if INCLUDE:
								if USEFILE:
									seq= SEQ[START:STOP]+seq
									aa= FILEseq[(ID,(START,STOP,STRAND))]+aa
								elif EXTEND:
									part= revcomp(SEQ[START:STOP])
									for i in range(0,len(part)-2,3):
										coding=part[i:i+3].upper()
										seq+= coding
										aa+=CODE.get(coding,'X')
									if len(part)%3!=0:
										warning+= ' CDS:NOT3x'
										if len(part)%3==2:
											aa+='2'
											seq+=part[-2:]
										else:	
											aa+='1'
											seq+=part[-1]					
								else:
									seq= SEQ[START:STOP]+seq			
							if DW>0 and not ((USEFILE and FILEseq[(ID,(START,STOP,STRAND))][-1]=='*') or (not USEFILE and CODE.get(revcomp(SEQ[START:START+3]),0)=='*')):
								if START-DW<0:
									warning+= ' DW:'+str(START)
									part= SEQ[:START]
								else:
									part= SEQ[START-DW:START]
								if EXTEND:
									s1=START+1
									part= revcomp(part)
									for i in range(0,len(part)-2,3):
										s1-=3
										coding=part[i:i+3].upper()
										seq+= coding
										aa+=CODE.get(coding,'X')
										if CODE.get(coding,0)=='*':
											break
								else:
									seq= part+seq
							else:
								s1=START+1
							if not EXTEND: #automatically excludes USEFILE
								seq= revcomp(seq)
						name= '>'+ID+' | '+str(START+1)+':'+str(STOP)+'['+STRAND+'] -> '+UgD
						if warning!='':
							warning= ' (WARNING:'+warning+')'
						if EXTEND:
							name+=' (ext. '+str(s1)+':'+str(s2)+')'
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
						#print name+'\t'+warning
						#print '    nt.\t'+seq
						#if EXTEND:
						#	print '    aa.\t'+aa
			if os.path.isfile(tmp):
				os.remove(tmp) #previously overwritten, needs to be removed after last iteration
			if BAD!=[] or BAD2!=[]:
				print
				print '*** UNMATCHED ID(s) WARNING ***'
				if BAD!=[]:
					for ID in BAD:
						print 'no match in nucleotide sequence file for ID:',ID
					print
				if BAD2!=[]:
					for pair in BAD2:
						ID=pair[0]
						start=str(pair[1][0]+1)
						stop= str(pair[1][1])
						strand=pair[1][2]
						print 'no match in (optional) protein file for ID:',ID,'\tat coordinates '+start+':'+stop+' ['+strand+']'
			print

