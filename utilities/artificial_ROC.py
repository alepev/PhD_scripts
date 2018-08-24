#! /usr/bin/python


cutoffL= range(101,1)
cutoffL.sort(reverse=True)
ROC={}
for cutoff in cutoffL:
	for pat in testL:
		if PROB[pat]>= cutoff:
			if patD[pat]==0:
				TP+=1
			else:
				FP+=1
		else:
			if patD[pat]==0:
				FN+=1
			else:
				TN+=1
	SENS,1_SPEC= TP/(TP+FN),TN/(TN+FP)
	ROC[cutoff]= (SENS,1_SPEC)
data_points= ROC.values()
data_points.sort()

for point in data_points:
	print point
