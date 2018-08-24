#!/usr/bin/python

import sys
ignore= True #default is to ignore inclusion threshold

if len(sys.argv) != 4:
	print '\nusage: <HMMER output> <fasta file> <score threshold>'
	print '\nscore considered is the one reported in column 2 of hmmsearch output file;'
	print 'first check the output to decided at which score you would like to cut the results'
	print 'all detected locus models are included in the results,'
	print 'however standard output will indicate BEST one and all those (in order) above score threshold'
	print '(decision was to include all model because score is just an indication,'
	print ' the most appropriate model should be chosen based on alignment with orthologous sequences)'
else:
	hmmer= sys.argv[1]
	fasta= sys.argv[2]
	S=float(sys.argv[3])
	out= '.'.join(hmmer.split('/')[-1].split('.')[:-1])+'.'+str(S)+'.fa'

	#load hmmsearch results

	HITS= {}
	model= {}
	MULTI=False
	for line in open(hmmer).readlines()[17:]:
		x= line.split()
		if len(x)>=9:
			score=float(x[1])
			ID=x[8]
			if score >= S:
				m= ID.split('|')[0].split('.')
				if len(m)>1 and len(m[-1])==1:
					m='.'.join(m[:-1])
					if model.get(m,0)!=0:
						model[m]= model[m]+[ID]
						MULTI=True
					else:
						model[m]=[ID]
				HITS[ID]=1
			else:
				break
		elif len(x)==4: #inclusion threshold
			if ignore:
				continue
			else:
				break
		else:
			break

	#parse fasta seqs
	f=open(out,'w')
	APPEND=False
	for line in open(fasta).readlines():
		if line[0]=='>':
			locus=line.split()[0][1:]
			if HITS.get(locus,0)!=0:
				HITS[locus]=2
				f.write(line)
				APPEND=True
			else:
				APPEND=False				
		elif APPEND:
			f.write(line)
	f.close()

	#output check (screen output)
	if 1 in HITS.values():
		print '***WARNING: no fasta sequence retrieved for the following hits:***'
		for locus in HITS.keys():
			if HITS[locus]!=2:
				print locus
	print

	#warning for multiple models
	if MULTI:
		print 'best model in case of multiple hits per locus:'
		for k in model.keys():
			if len(model[k])>1:
				print 'BEST:\t'+model[k][0]
				for m in model[k][1:]:
					print '\t'+m
		print
