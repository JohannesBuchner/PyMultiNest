import pymultinest
import math

def myprior(cube, ndim, nparams):
	#print "cube before", [cube[i] for i in range(ndim)]
	for i in range(ndim):
		cube[i] = cube[i] * 10 * math.pi
	#print "cube after", [cube[i] for i in range(ndim)]

def myloglike(cube, ndim, nparams):
	chi = 1.
	#print "cube", [cube[i] for i in range(ndim)], cube
	for i in range(ndim):
		chi *= math.cos(cube[i] / 2.)
	#print "returning", math.pow(2. + chi, 5)
	return math.pow(2. + chi, 5)

multinest.run(myloglike, myprior, 2, resume = False, verbose = True)






