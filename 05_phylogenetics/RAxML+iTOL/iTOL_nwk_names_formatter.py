#!/usr/bin/python

import sys

tree=	sys.argv[1]

T= open(tree).readlines()[0] #tree is in the 1st (and only) line of the file

S=T.replace("'",'')
T=S

S=''
i=0
while i<len(T):
	if T[i] not in [':','(',')',',',';','\n']:
		if len(T[i:].split(',')[0]) < len(T[i:].split(')')[0]):
			block= T[i:].split(',')[0]
		else:
			block= T[i:].split(')')[0]
		L=len(block)
		b1=block.split(':')[0]
		b2=block.split(':')[-1]
		if len(b1.split('|'))>1:
			if len(b1.split('|')[1].split('bZIP'))>1:
				b1= b1.split('|')[0]+'__'+b1.split('|')[1]
			elif b1[:2]=='gi':
				b1='|'.join(b1.split('|')[:4])
			else:
				b1=b1.split('|')[0]
		S=S+b1+':'+b2+T[i+L]
		i=i+L+1
	elif T[i]==':':
		if len(T[i:].split(',')[0]) < len(T[i:].split(')')[0]):
			block= T[i:].split(',')[0]
		else:
			block= T[i:].split(')')[0]
		S=S+block
		i=i+len(block)
	else:
		S+=T[i]		
		i+=1
print S	
