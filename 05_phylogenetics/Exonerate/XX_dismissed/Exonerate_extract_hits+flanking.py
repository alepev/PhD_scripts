#!/usr/bin/python

import sys
import subprocess as sub
sub.PIPE=1

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
		seq2+= compl.get(seq[i].upper(),'X')
	return seq2
	
#ABOUT EXONERATE INDEXING:
#in-between coordinate system, positions are counted between the symbols rather than on the symbols 
#ex. for the sequence "ACGT":
#     A C G T
#    0 1 2 3 4
#
#ABOUT OUTPUT INDEXING: 
#- always refer to [+] strand (i.e. original nt sequence)
#- residue indexing starts from 0 
#- limits: start inclusive, end exclusive
#ex. C...ATGCGT...TAG	
#    0   123456   789
#-> if ATGCGT is the aligned region, then limits are 1:7 (1 inclusive, 7 exclusive)

if len(sys.argv)<3 or len(sys.argv)>4:
	print '\nusage: <Exonerate output> <S> (0 for ALL HITS) <N+|-> (optional)'
	print '\n#mode:'
	print ' - score only:\n  extract AA/NT sequences for all hits above score\n  (S=0 for ALL HITS)'
	print ' - flanking region:\n  extract N dna flanking residues (upstream/downstream) for hits above score'
	print '\tex. 2000- == 5-prime / upstream'
	print '\tex. 2000+ == 3-prime / downstream'
	print
else:
	infile= sys.argv[1]
	S= sys.argv[2]
	if len(sys.argv)==4:
		F=sys.argv[3]
		FLANK= True
	else:
		FLANK=False

	if not S.isdigit() or int(S)<0 or (len(sys.argv)==4 and (not F[:-1].isdigit or F[-1] not in ['+','-'])):
		print '***ERROR - PLEASE CHECK YOUR INPUT!***'
	else:
		genomefile= open(infile).readlines()[0].split(' -t ')[1].split()[0] #genome file name
		S= int(S)
		QDB={}
		TDB={}
		BEST=[]
		uniqID=0
		RESULTS= '\n'.join(open(infile).readlines()).split('Query:')
		if len(RESULTS)>1:
			for q in RESULTS[1:]:
				qname= q.split()[0]
				#target sequence features
				tname= q.split('Target: ')[1].split('\n')[0].split(':[revcomp]')
				if len(tname)>1:
					strand='-'
				else:
					strand='+'
				tname=tname[0].split('len=')[0]
				score= int(q.split('Raw score: ')[1].split()[0])
				if score<S:
					continue
				x= q.split('Target range: ')[1].split(' -> ')
				start= int(x[0].split()[-1])
				stop= int(x[1].split()[0])
				if stop<start:
					tmp=stop
					stop=start
					start=tmp
				TDB[tname]=TDB.get(tname,[])+[(start,stop,strand,qname,score)]
				if not FLANK: #ONLY FOR AA/NT SEQUENCE EXTRACTION
					uniqID+=1
					BEST+= [(score,qname,tname,uniqID)]
					x= x[1].split('vulgar:')[0].split('\n')[1:-1]
					##target sequence
					aaseq=[]
					ntseq=[]
					for line in x:
						if len(line)>0:
							x1= line.lstrip()[0]
							if x1.isalpha() or x1 in ['*','-']:
								aaseq.append(line.lstrip().rstrip())
							elif x1.isdigit():
								ntseq+=[line.lstrip()]#.split()[3]])
					#nt sequence
					NTseq=''
					for i in range(len(ntseq)):
						if (i+1)%2 == 0: #only even lines are the ones to keep, the others contain query AA seq
							NTseq+=ntseq[i].split()[2]
					NTseq=NTseq.replace('-','') #remove gaps
					#aa sequence
					aaseq=''.join(aaseq).replace('-','') #remove gaps
					i=0
					c=0
					AAseq=''
					while i <(len(aaseq)-2):
						if aaseq[i].isupper():
							AAseq+=AA[aaseq[i:i+3]]
							i+=3
							c+=1
						elif aaseq[i]=='*':
							AAseq+='*'
							i+=3
						else:
							i+=1
					QDB[qname]=QDB.get(qname,[])+[(score,tname,start,stop,strand,c,AAseq,NTseq)]
		##FLANKING SEQUENCE EXTRACTION - OUTPUT
		if FLANK:
			N= int(F[:-1])
			if F[-1]=='-':
				UP=True
			else:
				UP=False
			#target sequences/contigs
			T= TDB.keys()
			for k in T:
				t= k.split()[0]
				if t[0]=='>':
					t=t[1:]
				x= TDB[k]
				TDB.pop(k)
				TDB[t]= x
			#retrieve selected target sequences/contigs
			SEQDB={}
			for line in open(genomefile).readlines():
				if len(line)>0 and line[0]=='>':
					name= line.split()[0][1:]
					if TDB.get(name,0)!=0:
						READ=True
					else:	
						READ=False
				elif READ:
					SEQDB[name]= SEQDB.get(name,'') + line.split()[0]
			for t in TDB.keys():
				if t not in SEQDB.keys():
					print 'not found ',t
			#output flanking sequence TDB[tname]=TDB.get(tname,[])+[(start,stop,strand,qname,score)]
			T= TDB.keys() #overwrites previous list
			T.sort()
			for t in T:
				sorder= TDB[t]
				sorder.sort()
				for hit in sorder:
					start= hit[0]
					stop= hit[1]
					strand = hit[2]
					warning= ''
					if strand=='+':
						if UP:
							if start-N<0:
								warning= ' (WARNING: extracted '+str(start)+' only)'
								seq= SEQDB[t][:start]
							else:
								seq= SEQDB[t][start-N:start]
						else:
							if stop+N>len(SEQDB[t]):
								warning= ' (WARNING: extracted '+str(len(SEQDB[t])-stop)+' only)'
								seq= SEQDB[t][stop:]
							else:
								seq= SEQDB[t][stop:stop+N]
					else:
						if UP:
							if stop+N>len(SEQDB[t]):
								warning= ' (WARNING: extracted '+str(len(SEQDB[t])-stop)+' only)'
								seq= SEQDB[t][stop:]
							else:
								seq= SEQDB[t][stop:stop+N]
						else:
							if start-N<0:
								warning= ' (WARNING: extracted '+str(start)+' only)'
								seq= SEQDB[t][:start]
							else:
								seq= SEQDB[t][start-N:start]
						seq= revcomp(seq)
					name= '>'+t+' | EX. vs '+hit[3]+' S'+str(hit[4])+' flank'+F+' ['+strand+'] '+str(start)+':'+str(stop)+warning
					print name
					print seq
			print

		##AA/NT SEQUENCE EXTRACTION - OUTPUT
		else:
			print
			print 'Exonerate output file:\t'+infile
			print 'reference genome file:\t'+genomefile
			print 
			print '****************************************************************'
			print
			#output - summary by score
			if len(BEST)>0:
				print '###SUMMARY - HITS BY SCORE###'
				print '#SCORE\tQUERY\tTARGET'
				BEST.sort(reverse=True)
				for b in BEST:
					print b[0],'\t',b[1],'\t\t',b[2]
			#output - summary by target
				print
				print '###SUMMARY - HITS BY TARGET###'
				print '#TARGET\tSTART\tSTOP\tSTRAND\tQUERY\tSCORE'
				T= TDB.keys()
				T.sort()
				for t in T:
					sorder= TDB[t]
					sorder.sort()
					for hit in sorder:
						print t+'\t'+str(hit[0])+'\t'+str(hit[1])+'\t'+hit[2]+'\t'+hit[3]+'\t'+str(hit[4])
			#output - details by query
				print 
				print '****************************************************************'
				print '###RESULTS DETAIL (BY QUERY)###'
				qs= QDB.keys()
				qs.sort()
				for q in qs:
					print
					print '#QUERY: ',q
					ss=QDB[q]
					ss.sort(reverse=True)
					for s in ss:
						score=s[0]
						tname=s[1]
						start=s[2]
						stop=s[3]
						strand=s[4]
						aalen=s[5]
						print
						print '>'+tname
						print s[-2]
						print s[-1]
						print '#SCORE:  ',score
						print '#LIMITS: ',str(start)+'-'+str(stop),'['+strand+']'
						print '#LENGTH:  AA',aalen,' DNA',stop-start
						if s[-2].count('M')==0:
							print '*WARNING: missing Start codon!\t\t\t(suggested flanking sequence extension)'
						elif s[-2][0]!='M':
							print '*WARNING: sequence not starting with Met!\t(correct Start codon might be inside or outside current sequence limits)'
						if s[-2].count('*')>0:
							print '*WARNING: Stop codon inside sequence!\t\t(might indicate a false positive hit)'
						print 
		







