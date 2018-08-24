#!/usr/bin/python

#recursive combinatorial function
#given a list of elements, creates a list of all possible groupings of elements (without repetition) 
#from length 2 (== pairs) to length equal to list length - 1

def CXl(PRE,l,i,C,R):
	if len(PRE)==C:
		R.append(PRE)
	else:
		for j in range (i+1,len(l)):
			CXl(PRE+[l[j]],l,j,C,R)

l= ['A','B','C','D','E','F']
C=2
R=[]

while C < len(l):
	for i in range(len(l)-C+1):
		PRE=[l[i]]
		CXl(PRE,l,i,C,R)
	C+=1

for r in R:
	print r
print len(R)
