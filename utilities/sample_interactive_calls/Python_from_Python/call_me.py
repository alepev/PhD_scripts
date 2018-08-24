#! /usr/bin/python

import sys
import time
import random
r= sys.argv[1]
time.sleep(random.choice([1,2,3,4,5]))
print 'started..',r
time.sleep(random.choice([1,2,3,4,5]))
print '..done!',r
