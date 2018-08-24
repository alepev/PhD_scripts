#!/usr/bin/python

MEME_COLORS= ['#00ffff','#0000ff','#ff0000','#ff00ff','#ffff00','#00ff00','#008080','#444444','#008000','#c0c0c0','#800080','#808000','#000080','#800000','#e2fffe']
Pvalue= 1e-10

#MOTIFS up to 15th:
#	- thick rectangle (RE) for those with p-value <= Pvalue, thin rectangle (GP) others
#MOTIFS > 15th and up to 30th:
#	- above colors are recycled
#	- ellipse (EL) for those with p-value <= Pvalue, diamond (DI) others

#motif syntax: shape|start|end|color|motif_name
#annotation syntax: sequence_name,length,motif1,motif2,..

#NOTE: sequence name MUST correspond to the one in nwk tree!!
#MEME may truncate the name if too long, correct it to allow matching

import sys,copy

if len(sys.argv)!= 2:
	print 'usage: <MEME.txt file>'
	print 'produces motif annotation to be associated to iTOL tree (MAX.30 MOTIFS!)'
else:
	MEME= sys.argv[1]
	out= '.'.join(MEME.split('.')[:-1])+'.ann'
	
	LM={}
	X= ''.join(open(MEME).readlines()).split('***\nMOTIF')[1:]
	for i in range(len(X)):
		w= int(X[i].split()[3])
		LM[i+1]=w
	X= ''.join(open(MEME).readlines()).split('SUMMARY OF MOTIFS')[1].split('\n')[8:]
	i=0
	while X[i][0]!='-':
		x=X[i].split()
		ID= x[0]
		BLOCKS= x[2]
		while BLOCKS[-1]=='\\':
			i+=1
			BLOCKS= BLOCKS[:-1]+X[i].split()[0]
		end=0
		ann=''
		for b in BLOCKS.split('_'):
			if b[0]!='[':
				end+= int(b)
			else:
				p= float(b.split('(')[1][:-2])
				N= int(b.split('(')[0][1:])
				if N<=len(MEME_COLORS):
					color= MEME_COLORS[N-1]
					if p<=Pvalue:
						shape= 'RE'
					else:
						shape= 'GP'
				else:
					color= MEME_COLORS[N-len(MEME_COLORS)-1] #colors recycled
					if p<=Pvalue:
						shape= 'EL'
					else:
						shape= 'DI'
				start=end+1
				end=end+LM[N]
				
				name= 'MOTIF_'+str(N)
				ann+= ','+'|'.join([shape,str(start),str(end),color,name])
		print ID+','+str(end)+ann
		i+=1
