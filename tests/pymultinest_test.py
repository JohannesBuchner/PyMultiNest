#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
import numpy
from numpy import pi, cos
from pymultinest.solve import solve
import json
import os
if not os.path.exists("chains"): os.mkdir("chains")

# probability function, taken from the eggbox problem.

def test_1():
	def myprior(cube):
		return cube * 10 * pi

	def myloglike(cube):
		chi = (cos(cube / 2.)).prod()
		return (2. + chi)**5

	# number of dimensions our problem has
	parameters = ["x", "y"]
	n_params = len(parameters)

	# run MultiNest
	result = solve(LogLikelihood=myloglike, Prior=myprior, 
		n_dims=n_params, outputfiles_basename="chains/3-")
	
	with open('%sparams.json' % "chains/3-", 'w') as f:
		json.dump(parameters, f, indent=2)

	print()
	print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
	print()
	print('parameter values:')
	for name, col in zip(parameters, result['samples'].transpose()):
		print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))

	assert numpy.isclose(result['logZ'], 236.0, atol=1), result['logZ']
	assert numpy.isclose(result['logZerr'], 0.1, atol=1), result['logZerr']

def test_2():
	"""testing n_params != n_dims case"""
	def myprior(cube):
		cube[-1] = numpy.sum(cube[:-1])
		return cube * 10 * pi

	def myloglike(cube):
		chi = (cos(cube[:-1] / 2.)).prod()
		return (2. + chi)**5

	# number of dimensions our problem has
	parameters = ["x", "y"]
	n_dims = len(parameters)
	n_params = n_dims + 1

	# run MultiNest
	result = solve(LogLikelihood=myloglike, Prior=myprior,
		n_dims=n_dims, n_params=n_params, outputfiles_basename="chains/5-")

	with open('%sparams.json' % "chains/5-", 'w') as f:
		json.dump(parameters, f, indent=2)

	parameters += ['x+y']
	print()
	print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
	print()
	print('parameter values:')
	for name, col in zip(parameters, result['samples'].transpose()):
		print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))

	assert numpy.isclose(result['logZ'], 239.1, atol=1), result['logZ']
	assert numpy.isclose(result['logZerr'], 0.1, atol=1), result['logZerr']



if __name__ == '__main__':
	test_1()
	test_2()

