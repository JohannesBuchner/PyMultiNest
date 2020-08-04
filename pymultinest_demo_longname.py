from __future__ import absolute_import, unicode_literals, print_function
import pymultinest
import math, os
try: os.mkdir('chains')
except OSError: pass

# Based in pymultinest_demo_minimal, this tests MultiNest when running with
# a long file name (> than 100 characters)

# our probability functions
# Taken from the eggbox problem.

# Take a look at pymultinest_demo.py for more pythonic convenience

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

outputfiles_basename="chains/a_really_really_really_really_really_really_really_really_really_really_really_really_really_really_really_really_long_name-"

# run MultiNest
pymultinest.run(myloglike, myprior, n_params, 
	resume = True, verbose = True,
    outputfiles_basename = outputfiles_basename)

# for making marginal / corner plots and reading the output files --> see full demo


