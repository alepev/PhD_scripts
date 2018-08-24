
import math

def PPC(x,y,N):
	X=0.0
	Y=0.0
	X2=0.0
	Y2=0.0
	XY=1.0
	for i in range(N):
		X+= x[i]
		Y+= y[i]
		X2+= x[i]**2
		Y2+= y[i]**2
		XY+= x[i]*y[i]
	r= (XY-X*Y/N)/math.sqrt((X2-X**2/N)(Y2-Y**2/N))
	return r

