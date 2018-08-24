#!/usr/bin/python

#v2.0 keeps full sequence name, while previously cut in case of spacing characters

import sys

if len(sys.argv)<2 or len(sys.argv)>3:
	print
	print 'usage: <fasta file> <r|s> (optional, to create new file with removed duplicates / check for subsequences)'
	print
	print '	- DEFAULT: checks presence of 100% identical sequences (handles presence of final * and mixed lower/uppercase)'
	print '	- R option: removes duplicates (1st sequence in alphabetical order is the one kept)'
	print '	- S option: additionally checks presence of sequences with 100% identical aligned regions,'
	print '	            allowing for insertion/deletions but not mismatches'
	print '	            (NOTE: USE ONLY WITH ALIGNED SEQUENCES, otherwise it does`t make sense!)' 
	print
else:
	fasta= sys.argv[1]
	remove=False
	subseq=False
	if len(sys.argv)==3:
		if sys.argv[2].upper()=='R':
			remove=True
			outfile= '.'.join(fasta.split('.')[:-1])+'.clean.fa'
			logfile= '.'.join(fasta.split('.')[:-1])+'.clean.log'
		elif sys.argv[2].upper()=='S':
			subseq=True
			subsDB={}
		else:
			print 
			print '***WARNING: unknown argument "'+sys.argv[2]+'" (ignored)***'

	#load fasta 
	seqDB={}
	ID= open(fasta).readlines()[0][:-1]
	seq=''
	for line in open(fasta).readlines()[1:]:
		if len(line.split())>0:
			if line[0]=='>':
				seq=seq.upper()
				if seq[-1]=='*':
					seq=seq[:-1]
				seqDB[ID]=seq	
				ID=line[:-1]
				seq=''
			else:
				seq+= line.split()[0]
	seq=seq.upper()
	if seq[-1]=='*':
		seq=seq[:-1]
	seqDB[ID]=seq	
	
	#check for duplicates
	sameDB={}
	dupli={}
	DB=seqDB.keys()
	DB.sort()
	for s in range(len(DB)-1):
		s1= DB[s]
		if dupli.get(s1,0)==0:
			for s2 in DB[s+1:]:
				if dupli.get(s2,0)==0 and len(seqDB[s1])==len(seqDB[s2]):
					seq1=seqDB[s1]
					seq2=seqDB[s2]
					EQUAL=True
					if subseq:
						gap1,gap2=0,0
					for i in range(len(seq1)):
						if seq1[i]!=seq2[i]:
							if subseq:
								if seq1[i]=='-':
									gap1+=1
									continue
								elif seq2[i]=='-':
									gap2+=1
									continue
							EQUAL=False
							break
					if EQUAL:
						if not subseq or (subseq and gap1+gap2==0):
							sameDB[s1]=sameDB.get(s1,[])+[s2]
							dupli[s2]=1
						else:
							subsDB[s1]= subsDB.get(s1,[])+[(gap1,gap2,s2)]
	
	#output info
	print
	if len(sameDB)>0:
		if remove:
			f=open(logfile,'w')
			f.write('***duplicate sequences*** (only 1st kept)')
		print '***found duplicate sequences***'
		DB=sameDB.keys()
		DB.sort()
		for k in DB:
			if remove:
				f.write('\n\n'+k)
			print
			print k
			for k2 in sameDB[k]:
				if remove:
					f.write('\nduplicate:\t'+k2)
				print 'duplicate:\t'+k2
		print
		print '**************************'
		if remove:
			f.close()
			DB=seqDB.keys()
			DB.sort()
			f=open(outfile,'w')
			for k in DB:
				if dupli.get(k,0)==0:
					f.write(k+'\n'+seqDB[k]+'\n')
			f.close()
			print 'NON-REDUNDANT FILE CREATED: (see log file for details)\n'+outfile
		else:
			print '(add the optional argument "r" to generate a non-redundant version of the file)'
	else:
		print '(no duplicates detected)'	
	print
	if subseq:
		if len(subsDB)>0:
			print '***found identical sequences with indels***'
			DB=subsDB.keys()
			DB.sort()
			for k in DB:
				print 
				print k
				for k2 in subsDB[k]:
					print 'alternative sequence:\t'+k2[2]
					print 'gaps in main:\t'+str(k2[0])
					print 'gaps in alternative:\t'+str(k2[1])
		else:
			print '(no indels detected)'	
		print
