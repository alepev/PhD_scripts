#!/usr/bin/python

gff3='TAIR9_GFF3_genes.gff'

X={}
for line in open(gff3,'rU').readlines()[1:]:
	x=line.split()[2]
	X[x]=X.get(x,0)+1
K=X.keys()
K.sort()
for k in K:
	print X[k],k

'''
ex. of results with TAIR9_GFF3_genes.gff:
179783 CDS
6 chromosome
196957 exon
31257 five_prime_UTR
28691 gene
37318 mRNA
176 miRNA
428 ncRNA
33410 protein
926 pseudogene
1263 pseudogenic_exon
930 pseudogenic_transcript
15 rRNA
13 snRNA
71 snoRNA
689 tRNA
28076 three_prime_UTR
3901 transposable_element_gene
+
31189 transposable_element
34856 transposon_fragment
'''