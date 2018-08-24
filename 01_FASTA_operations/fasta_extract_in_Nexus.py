#! /usr/bin/python

fasta='ALL+ref_S1fa'
nexus='S1_bZIP1+53.txt'
parsed='.'.join(nexus.split('.')[:-1])+'.fa'

PARSE={}

for line in open(nexus).readlines()[4:]:
	if line[0]==';':
		break
	else:
		name=line.split()[1].split(',')[0].replace('=',':') #removes initial number and final comma
		PARSE[name]=1
		
f= open(parsed,'w')
for line in open(fasta):
	x=line.split()[0]
	if line[0]=='>':
		name=x[1:]
	 	if PARSE.get(name,0)!=0:
			f.write(line)
			READ=True
			PARSE[name]=2
		else:
			READ=False
	elif READ:
		f.write(line)
f.close()
		
if any(map(lambda x: PARSE[x]==1, PARSE.keys())):
	print 'following sequence(s) not found!'
	for name in PARSE.keys():
		if PARSE[name]==1:
			print name
		
	