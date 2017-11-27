from __future__ import absolute_import, unicode_literals, print_function
import pymultinest
import math, os
if not os.path.exists("chains"): os.mkdir("chains")

# This "old" api may be a little bit faster and suitable for likelihoods
# written in C (via ctypes).

# our probability functions
# Taken from the eggbox problem.

# Take a look at pymultinest_demo.py for more pythonic convenience
def test():

	def myprior(cube, ndim, nparams):
		for i in range(ndim):
			cube[i] = cube[i] * 10 * math.pi

	def myloglike(cube, ndim, nparams):
		chi = 1.
		for i in range(ndim):
			chi *= math.cos(cube[i] / 2.)
		return math.pow(2. + chi, 5)

	# number of dimensions our problem has
	parameters = ["x", "y"]
	n_params = len(parameters)

	# run MultiNest
	pymultinest.run(myloglike, myprior, n_params, 
		resume = True, verbose = True)

	
