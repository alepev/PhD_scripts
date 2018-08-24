#!/usr/bin/python

import sys

def any(iterable):
	for element in iterable:
		if element:
			return True
	return False

if len(sys.argv)!= 2:
	print 'usage: <motif list>'
else:
	sitesfile= sys.argv[1]

	print '#original\t#IUPAC\t#rev.complementary(IUPAC)'

	#consensus code
	IUPAC={}					
	IUPAC['N']=['A','C','G','T']
	IUPAC['R']=['A','G']
	IUPAC['Y']=['C','T']
	IUPAC['S']=['C','G']
	IUPAC['W']=['A','T']
	IUPAC['K']=['G','T']
	IUPAC['M']=['A','C']
	IUPAC['B']=['C','G','T']
	IUPAC['D']=['A','G','T']
	IUPAC['H']=['A','C','T']
	IUPAC['V']=['A','C','G']
	inv_IUPAC={}
	for key in IUPAC.keys():
		inv_IUPAC[str(IUPAC[key])]=key
	#complementary nt
	compl={'A':'T','C':'G','G':'C','T':'A'} 
	compl['N']='N'
	compl['R']='Y'
	compl['Y']='R'
	compl['S']='S'
	compl['W']='W'
	compl['K']='M'
	compl['M']='K'
	compl['B']='V'
	compl['V']='B'
	compl['D']='H'
	compl['H']='D'

	sitelist=[]
	reverse={}
	for line in open(sitesfile).readlines():
		if line[0] in compl.keys() or line[0]=='(':
			x=line.split()
			if x[0].count('(')==0:
				d=x[0]
			else:
				d=''
				i=0
				while i < len(x[0]):
					if x[0][i]!='(':
						d+=x[0][i]
						i+=1
					else:
						part=x[0][i+1:].split(')')
						comb=part[0].split('/')
						comb.sort()
						d+=inv_IUPAC[str(comb)][0]
						i=i+len(part[0])+2
			if d not in sitelist and d not in reverse.values():
				sitelist.append(d)
				r=map(lambda x: compl[x], d)
				r.reverse()
				r= ''.join(r)
				reverse[d]=r
				print x[0],'\t',d,'\t',r
