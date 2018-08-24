#!/usr/bin/python

import sys

topfile= 'TopALL.mapped'

#extracts and visualize loci and their expression values from a list of UP-/DW-regulated genes
#(NOTE: change it to work also with topTables! should be able to extract expr. for one locus - all conditions if none specified)
if len(sys.argv)!=2:
	print 'usage: <gene_list> (1 column listing loci, no headers) --> SCREEN OUTPUT'
else:	
	listfile=str(sys.argv[1])
	
	L=[]
	for line in open(listfile).readlines():
		x=line.split()
		L.append(x[0])
	print open(topfile).readlines()[:2]
	for line in open(topfile).readlines()[2:]:
		if line.split()[0] in L:
			print line,
	#NOTE: comma used since printed lines already include \n character
	
