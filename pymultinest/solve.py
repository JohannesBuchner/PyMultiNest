"""
A pythonic interface to MultiNest
"""

from __future__ import absolute_import, unicode_literals, print_function

from .run import run
from .analyse import Analyzer
import numpy
import tempfile
import shutil

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
def solve(LogLikelihood, Prior, n_dims, **kwargs):
	kwargs['n_dims'] = n_dims
	files_temporary = False
	if 'outputfiles_basename' not in kwargs:
		files_temporary = True
		tempdir = tempfile.mkdtemp('pymultinest')
		kwargs['outputfiles_basename'] = tempdir + '/'
	outputfiles_basename = kwargs['outputfiles_basename']
	def SafePrior(cube, ndim, nparams):
		try:
			n_params = kwargs.get('n_params', n_dims)
			a = numpy.array([cube[i] for i in range(n_params)])
			b = Prior(a)
			for i in range(n_params):
				cube[i] = b[i]
		except Exception as e:
			import sys
			sys.stderr.write('ERROR in prior: %s\n' % e)
			sys.exit(1)
	
	def SafeLoglikelihood(cube, ndim, nparams, lnew):
		try:
			a = numpy.array([cube[i] for i in range(n_dims)])
			l = float(LogLikelihood(a))
			if not numpy.isfinite(l):
				import sys
				sys.stderr.write('WARNING: loglikelihood not finite: %f\n' % (l))
				sys.stderr.write('         for parameters: %s\n' % a)
				sys.stderr.write('         returned very low value instead\n')
				return -1e100
			return l
		except Exception as e:
			import sys
			sys.stderr.write('ERROR in loglikelihood: %s\n' % e)
			sys.exit(1)
	
	kwargs['LogLikelihood'] = SafeLoglikelihood
	kwargs['Prior'] = SafePrior
	run(**kwargs)
	
	analyzer = Analyzer(n_dims, outputfiles_basename = outputfiles_basename)
	stats = analyzer.get_stats()
	samples = analyzer.get_equal_weighted_posterior()[:,:-1]
	
	if files_temporary:
		shutil.rmtree(tempdir, ignore_errors=True)
	
	return dict(logZ=stats['nested sampling global log-evidence'],
		logZerr=stats['nested sampling global log-evidence error'],
		samples = samples,
		)

class Solver(object):
	def __init__(self, **kwargs):
		self.n_dims = kwargs['n_dims']
		if 'outputfiles_basename' in kwargs:
			self.outputfiles_basename = kwargs['outputfiles_basename']
		else:
			self.outputfiles_basename = '(temporary directory)'
		results = solve(self.LogLikelihood, self.Prior, **kwargs)
		for k, v in results.items():
			setattr(self, k, v)
	def __str__(self):
		s = 'Model in "%s"' % self.outputfiles_basename
		s += ' (%d dimensions)\n' % self.n_dims
		s += 'Evidence ln Z = %.1f +- %.1f\n' % (self.logZ, self.logZerr)
		s += 'Parameter values:\n'
		for i, col in enumerate(self.samples.transpose()):
			s += '   Parameter %2d : %.3f +- %.3f\n' % (i + 1, col.mean(), col.std())
		return s.rstrip()
	def __repr__(self):
		return 'Solver(n_dims=%d, outputfiles_basename="%s")' % (self.n_dims, self.outputfiles_basename)


