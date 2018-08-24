#! /usr/bin/python

AA= {}
AA['Ala']='A'
AA['Arg']='R'
AA['Asn']='N'
AA['Asp']='D'
AA['Cys']='C'
AA['Glu']='E'
AA['Gln']='Q'
AA['Gly']='G'
AA['His']='H'
AA['Ile']='I'
AA['Leu']='L'
AA['Lys']='K'
AA['Met']='M'
AA['Phe']='F'
AA['Pro']='P'
AA['Ser']='S'
AA['Thr']='T'
AA['Trp']='W'
AA['Tyr']='Y'
AA['Val']='V'

compl={'A':'T','C':'G','G':'C','T':'A'} 
def revcomp(seq):
	seq2=''
	for i in range(len(seq)-1,-1,-1):
		seq2+= compl.get(seq[i].upper(),'N')
	return seq2

import sys
import os

if len(sys.argv)!=4:
	print
	print 'usage: <sequence_file> <limits_file> <5\'+seq+3\'>'
	print
	print '5\' + seq + 3\' extracts upstream + within-limits + downstream sequence'
	print '(5\', 3\' -> integer; seq -> 1/0 to include/leave out sequence within limits)'
	print
	print '***LIMITS FILE SPECIFICATIONS***'
	print ' - lines starting with "#" ignored'
	print ' - limits should be provided as seqNAME,START,STOP,(STRAND) [string,integer,integer,(+/-)];' 
	print '    if STRAND not given, START>STOP interpreted as reverse filament (STRAND = -)'
	print '    if given, START<STOP works with both +/-, however START>STOP will generate an error with STRAND = +'
	print ' - limits should be INCLUSIVE OF BOTH START/END residues, and indexing STARTS WITH 1'
	print '    during extraction limits are converted to allow string operations'
	print '    ex. 4-residue sequence starting from 1: sequence 1:4 -> corresponding to string 0:3 (in string notation: string[0:4])'
	print ' - seqNAME starting with ">" tolerated'
	print ' - seqNAME with additional descriptors tolerated if following "|", ex. ATXXXX|bZIPX (bZIPX ignored during ID matching)'
	print
else:
	#file names
	genfile= sys.argv[1]
	annfile= sys.argv[2]
	UgD= sys.argv[3]
	UP= abs(int(UgD.split('+')[0]))
	DW= abs(int(UgD.split('+')[2]))
	INCLUDE=False
	if UgD.split('+')[1]=='1':
		INCLUDE=True

	outfile= '.'.join(annfile.split('.')[:-1])+'_'+UgD+'.fa'

	#load limits
	DB={}
	ERROR=False
	for line in open(annfile).readlines():
		x=line.split('\t')
		if len(x)>=3 and x[0][0]!='#':
			ID= x[0].split()[0].split('|')[0]
			if ID[0]=='>':
				ID=ID[1:]
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
			DB[ID]=DB.get(ID,[])+[(START-1,STOP,STRAND)]
	
	if not ERROR:	
		#load fasta		
		SEQ={}	
		for line in open(genfile).readlines():
			if len(line.split())!=0:
				if line[0]=='>':
					ID= line[1:].split()[0].split('|')[0]
					if DB.get(ID,0)!=0:
						READ=True
					else:
						READ=False
				elif READ==True:
						SEQ[ID]=SEQ.get(ID,'')+line.split()[0]
		#estract selected sequence coordinates
		IDlist= DB.keys()
		IDlist.sort()
		BAD=[]
		f=open(outfile,'w')
		for ID in IDlist:
			if SEQ.get(ID,0)==0:
				BAD.append(ID)
			else:
				coord= DB[ID]
				coord.sort()
				for xy in coord:
					START= xy[0]
					STOP= xy[1]
					STRAND= xy[2]
					warning= ''
					seq=''
					print START+1,STOP,STRAND
					if STRAND=='+':
						if UP>0:
							if START-UP<0:
								warning+= 'UP:'+str(START)
								seq+= SEQ[ID][:START]
							else:
								seq+= SEQ[ID][START-UP:START]
						if INCLUDE:
							seq+= SEQ[ID][START:STOP]
						if DW>0:
							if STOP+DW>len(SEQ[ID]):
								warning+= 'DW:'+str(len(SEQ[ID])-STOP)
								seq+= SEQ[ID][STOP:]
							else:
								seq+= SEQ[ID][STOP:STOP+DW]
					else: #reverse strand [-]
						if DW>0:
							if START-DW<0:
								warning+= 'DW:'+str(START)
								seq+= SEQ[ID][:START]
							else:
								seq+= SEQ[ID][START-DW:START]
						if INCLUDE:
							seq+= SEQ[ID][START:STOP]
						if UP>0:
							if STOP+UP>len(SEQ[ID]):
								warning+= 'UP:'+str(len(SEQ[ID])-STOP)
								seq+= SEQ[ID][STOP:]
							else:
								seq+= SEQ[ID][STOP:STOP+UP]
						seq= revcomp(seq)		
					name= '>'+ID+' | '+str(START+1)+':'+str(STOP)+'['+STRAND+'] '+UgD
					if warning!='':
						name+= ' (WARNING: '+warning+')'
					f.write(name+'\n')
					f.write(seq+'\n')
					print name
					#print seq
		f.close()
		for ID in BAD:
			print 'no sequence retrieved for',ID
