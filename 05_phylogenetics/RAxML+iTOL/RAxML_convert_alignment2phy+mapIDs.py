#!/usr/bin/python

#v2.0: added replacement of commas (also cause an error in RAxML)

import sys
sortalign= False

align= sys.argv[1] #alignment in fasta/muscle format
phyout= '.'.join(align.split('.')[:-1])+'.phy'
logfile= '.'.join(align.split('.')[:-1])+'.phy.log'

SEQS={}
S=[] #keeps sequence order
MAP={}
multi= {}
for line in open(align).readlines():
	if line[0]=='>':
		name= line[1:-1]
		x=name.split()[0].split('/')[0]	#removes descriptions after spaces (if any) and sequence limits
		x=x.replace('[','_')
		x=x.replace(']','')
		x=x.replace(':','=')
		x=x.replace(';','-')
		x=x.replace(',','-')
		x=x.replace('(','_')
		x=x.replace(')','')
		S.append(name)
		multi[x]= multi.get(x,0)+1
		if multi[x]>1:
			MAP[name]=x+'_#'+str(multi[x])
		else:
			MAP[name]=x
	else:
		SEQS[name]= SEQS.get(name,'')+line[:-1] #removes '\n'

if sortalign:
	S.sort()

Nseq= len(SEQS.keys())
Lseq= len(SEQS[SEQS.keys()[0]])

MULTI=False
for v in multi.values():
	if v>0:
		MULTI=True
		break

f=open(phyout,'w')
f.write(str(Nseq)+' '+str(Lseq))
if MULTI:
	f2=open(logfile,'w')
for locus in S:
	f.write('\n'+MAP[locus]+' '+SEQS[locus])
	if MULTI:
		f2.write(MAP[locus]+'\t'+locus+'\n')
f.close()
if MULTI:
	f2.close()


