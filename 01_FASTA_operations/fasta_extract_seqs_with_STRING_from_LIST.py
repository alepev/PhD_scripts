#!/usr/bin/python

import sys
import subprocess as sub
sub.PIPE=1

if len(sys.argv)!=3:
	print 'usage: <fasta file> <string_list>'
	print 'extracts fasta sequence for which at least one of the strings appears in sequence name'
	print 'NOTE: '
	print '- only column1 of list file used'
	print '- strings should not contain spaces, otherwise truncated'
else:
	fasta= sys.argv[1]
	Sfile=sys.argv[2]
	
	out= '.'.join(fasta.split('/')[-1].split('.')[:-1])+'.'+Sfile.split('/')[-1].split('.')[0]+'_LIST.fa'
	
	#get STRING list
	
	print '\n***sequences list***'
	Slist=[]
	for line in open(Sfile,'rU').readlines():
		print line,
		S=line.split()[0]
		#S=S.split('|')[0] #removes description to ease matching
		S=S.split('/')[0] #removes sequence limits added by some programs, ex. Jalview
		if S!='':
			Slist.append(S)

	#call shell to output 

	print '\n***matches to fasta file***'
	OK={}
	for S in Slist:
		print '\nquery:\t',S,
		x=sub.Popen('grep "'+S+'" '+fasta, shell=True, stdout=1)
		headers= x.communicate()[0].split('\n')[:-1]
		if len(headers)==0:
			print '\t***NOT FOUND!***'
		else:
			print
			for line in headers:
				OK[line]=1
				print 'match:\t'+line
			print
	
	#parse fasta seqs

	f=open(out,'w')
	APPEND=True
	for line in open(fasta).readlines():
		if line[0]=='>':
			if OK.get(line[:-1])==1:
				f.write(line)
				APPEND=True
			else:
				APPEND=False				
		elif APPEND:
			f.write(line)
	f.close()

	print '\n***parsed sequences output file***'
	print out
	print

	
