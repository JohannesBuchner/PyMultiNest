#!/usr/bin/env python
from pymultinest.solve import Solver
from numpy import pi, cos

class EggboxProblem(Solver):
	def Prior(self, cube):
		return cube * 10 * pi
	def LogLikelihood(self, cube):
		chi = (cos(cube / 2.)).prod()
		return (2. + chi)**5

solution = EggboxProblem(n_dims = 2)

print(solution)
""" prints:
Model in "chains/4-" (2 dimensions)
dimensionality = 2
Evidence ln Z = 236.0 +- 0.1
Parameter values:
   Parameter 0 : 15.808 +- 7.874
   Parameter 1 : 14.840 +- 9.237
"""
# use: solution.logZ, solution.logZerr, solution.samples

