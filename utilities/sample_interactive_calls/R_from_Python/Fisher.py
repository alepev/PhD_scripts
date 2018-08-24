#!/usr/bin/python

#computes Fisher exact test calling R

import subprocess as sub
sub.PIPE=1

N=2000	#tot. genes in gene universe
K=200	#tot genes with annotation GOx
n=230	#D.E. genes
k=52	#D.E. genes with annotation GOx
x=sub.Popen('Rscript Fisher.R '+str(N)+' '+str(K)+' '+str(n)+' '+str(k), shell=True, stdout=1)
Fisher_p-value=x.communicate()[0]

print Fisher_p-value



