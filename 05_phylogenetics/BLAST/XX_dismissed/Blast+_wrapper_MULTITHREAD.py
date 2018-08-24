#!/usr/bin/python

import sys
import subprocess as sub
sub.PIPE=1

#MULTITHREAD OPTION DEACTIVATED BECAUSE OF TOO LONG TIME FOR DB GENERATION (MULTITHREADING SUPPORTED IN -DB BUT NOT IN -SUBJECT OPTION)
#NOTE: BLAST RUNS OUT OF MEMORY ON DESKTOP, USE MUTANT INSTEAD! 
Ncores= 1#int(sub.Popen('grep processor /proc/cpuinfo | wc -l', shell=True, stdout=1).communicate()[0])-1
if Ncores>1:
	MULTI=True
else:
	MULTI=False

BLAST= 'tblastn'#'blastp'#
E= 	0.01 #max. e-value - in theory 0.001, but if I put that I always get with lower than 0.001 only (seems a non-inclusive threshold)

#cfr with 
#time blast -p tblastn -i bZIP1_5UTR_peptide24.fa -j  picea_abies.master-rna-scaff.nov2012.fa -N T -e 0.01 -v 10 -a 12

if len(sys.argv)!=4:
	print '\nusage: <query> <target> <N>'
	print '\nDEFAULT SETTINGS:'
	print '- BLAST TYPE:\t'+BLAST
	print '- max e-value:\t'+str(E)
	print
else:
	query= sys.argv[1]
	subject= sys.argv[2]
	N= int(sys.argv[3])
	out= '.'.join(subject.split('/')[-1].split('.')[:-1])+'__'+'.'.join(query.split('/')[-1].split('.')[:-1])+'.out'
	print
	print '#GENERAL SETTINGS'
	print '- BLAST TYPE:\t',BLAST
	print '- max N hits:\t',N
	print '- max e-value:\t',E
	print

	#RUN BLAST
	if MULTI:
		db= '.'.join(subject.split('.')[:-1])+'.tmpdb'
		#convert subject to db - necessary for use of multithread option
		if BLAST=='tblastp':
			mkdb= 'makeblastdb -in '+subject+' -input_type fasta -dbtype prot -out '+db+''
		else: #blastn,tblastn,tblastx
			mkdb= 'makeblastdb -in '+subject+' -input_type fasta -dbtype nucl -out '+db+''
		sub.Popen(mkdb,shell=True, stdout=1).communicate()[0]
		#blast run
		blast= ''+BLAST+' -db '+db+' -query '+query+' -max_target_seqs '+str(N)+' -evalue '+str(E)+' -num_threads '+str(Ncores)+' -outfmt "6 qseqid sseqid score evalue sstart send sseq"' 
	else:
		blast= ''+BLAST+' -subject '+subject+' -query '+query+' -max_target_seqs '+str(N)+' -evalue '+str(E)+' -outfmt "6 qseqid sseqid score evalue sstart send sseq"' 
	#command debug - check how the command passed to the shell looks like
	#hits= sub.Popen('echo "'+blast+'"',shell=True, stdout=1).communicate()[0] 
	hits= sub.Popen(blast,shell=True, stdout=1).communicate()[0]
	print '#QUERY\tSUBJECT\tSCORE\tE-VALUE\tSUBJ.START\tSUBJ.END\tSUBJ.SEQ(AA)'
	print hits
	hits=hits.split('\n')
	if MULTI:
		#remove db
		rm= sub.Popen('rm '+db+'*',shell=True)

	#FILTER HITS
	HITS= {}
	for h in hits:
		break
		x= h.split('\t')
		if len(x)!=3:	
			continue	
		ID=x[0]
		score= int(x[1])
		if x[2][0]=='e':
			evalue= float('1'+x[2])
		else:
			evalue= float(x[2])
		if score >= S: #e-value filtering already in blast, not necessary here
			if HITS.get(ID,0)==0 or (MIN and (HITS[ID][1]<evalue or (HITS[ID][1]==evalue and HITS[ID][0]>score))) or (not MIN and (HITS[ID][1]>evalue or (HITS[ID][1]==evalue and HITS[ID][0]<score))): 
				HITS[ID]=(score,evalue) 
	hits= None #free memory

	if len(HITS.keys())==0:
		print '***NO SIGNIFICANT HITs!***'

	else:
		H= HITS.keys()
		H.sort()
		f=open(out,'w')
		for h in H:
			if len(h.split('|'))==1:
				ID= h.split()[0]
			else:
				ID= '|'.join(h.split('|')[:-1])
			values = '\t%d\t%.2e' % (HITS[h][0], HITS[h][1])
			f.write(ID+values+'\n')
		f.close() 

		print '\nfiltered hits are in',out
	print
