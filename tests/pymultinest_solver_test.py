#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
from pymultinest.solve import Solver
import numpy
from numpy import pi, cos

class EggboxProblem(Solver):
	def Prior(self, cube):
		return cube * 10 * pi
	def LogLikelihood(self, cube):
		chi = (cos(cube / 2.)).prod()
		return (2. + chi)**5

def test_solver():
	solution = EggboxProblem(n_dims = 2)

	print(solution)
	assert solution.n_dims == 2, solution.n_dims
	assert numpy.isclose(solution.logZ, 236.0, atol=1), solution.logZ
	assert numpy.isclose(solution.logZerr, 0.1, atol=1), solution.logZerr
	""" prints:
	Model in "chains/4-" (2 dimensions)
	dimensionality = 2
	Evidence ln Z = 236.0 +- 0.1
	Parameter values:
	   Parameter 0 : 15.808 +- 7.874
	   Parameter 1 : 14.840 +- 9.237
	"""
	# use: solution.logZ, solution.logZerr, solution.samples

if __name__ == '__main__':
	test_solver()
