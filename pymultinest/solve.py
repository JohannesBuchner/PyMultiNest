"""
A pythonic interface to MultiNest
"""

from __future__ import absolute_import, unicode_literals, print_function

from .run import run
from .analyse import Analyzer
import numpy

"""
A pythonic interface to MultiNest. The arguments are the same as in 
*run* (see there), except that the functions are defined differently for 
convenience.

@param Prior:
	Takes a numpy array and returns the transformed array.
	Example:
	def myprior(cube):
		return cube * 20 - 10

@param Loglikelihood:
	takes a numpy array of the transformed parameters,
	and should return the loglikelihood


@param nparams:
	dimensionality of the problem

"""
def solve(LogLikelihood, Prior, **kwargs):
	n_dims = kwargs['n_dims']
	outputfiles_basename = kwargs['outputfiles_basename']
	def SafePrior(cube, ndim, nparams):
		try:
			a = numpy.array([cube[i] for i in range(ndim)])
			b = Prior(a)
			for i in range(ndim):
				cube[i] = b[i]
		except Exception as e:
			import sys
			sys.stderr.write('ERROR in prior: %s' % e)
			sys.exit(1)
	
	def SafeLoglikelihood(cube, ndim, nparams, lnew):
		try:
			a = numpy.array([cube[i] for i in range(ndim)])
			l = float(LogLikelihood(a))
			if not numpy.isfinite(l):
				import sys
				sys.stderr.write('WARNING: loglikelihood not finite: %f' % (l))
				sys.stderr.write('         for parameters: %s' % a)
				sys.stderr.write('         returned very low value instead')
				sys.exit(1)
				return -1e100
			return l
		except Exception as e:
			import sys
			sys.stderr.write('ERROR in loglikelihood: %s' % e)
			sys.exit(1)
	
	kwargs['LogLikelihood'] = SafeLoglikelihood
	kwargs['Prior'] = SafePrior
	run(**kwargs)
	
	analyzer = Analyzer(n_dims, outputfiles_basename = outputfiles_basename)
	stats = analyzer.get_stats()
	
	return dict(logZ=stats['nested sampling global log-evidence'],
		logZerr=stats['nested sampling global log-evidence error'],
		samples = analyzer.get_equal_weighted_posterior()[:,:-1],
		)

	

