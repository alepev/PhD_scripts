#! /usr/bin/python

import subprocess as sub
import os
import time

Dpopen= {}
for i in range(5):
	Dpopen[i]= sub.Popen(["./call_me.py", str(i+1)])
while not all(map(lambda x: Dpopen[x].poll()==0,Dpopen.keys())): #WAITS FOR ALL THE PROCESSES TO BE CONCLUDED (checks every 2 secs)
	time.sleep(2)
print '***done***'
