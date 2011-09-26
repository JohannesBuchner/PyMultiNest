
from ctypes import cdll
lib = cdll.LoadLibrary('./libcnest.so')

from ctypes import *

def run(LogLikelihood,
	Prior,
	n_dims, 
	n_params = None, 
	n_clustering_params = None, wrapped_params = None, 
	multimodal = True, const_efficiency_mode = False, n_live_points = 1000,
	evidence_tolerance = 0.5, sampling_efficiency = 0.8, 
	n_iter_before_update = 100, null_log_evidence = -1e90,
	max_modes = 100,
	outputfiles_basename = "chains/1-", seed = -1, verbose = False,
	resume = True, context = 0):
	"""
	Run MultiNest
	
	@param n_params: 
		Total no. of parameters, should be equal to ndims in most cases 
		but if you need to store some additional
		parameters with the actual parameters then you need to pass 
		them through the likelihood routine.

	@param sampling_efficiency:
		defines the sampling efficiency. 0.8 and 0.3 are recommended 
		for parameter estimation & evidence evalutation
		respectively.


	@param evidence_tolerance:
		A value of 0.5 should give good enough accuracy.

	@param n_clustering_params:
		If mmodal is T, MultiNest will attempt to separate out the 
		modes. Mode separation is done through a clustering
		algorithm. Mode separation can be done on all the parameters 
		(in which case nCdims should be set to ndims) & it
		can also be done on a subset of parameters (in which case 
		nCdims < ndims) which might be advantageous as
		clustering is less accurate as the dimensionality increases. 
		If nCdims < ndims then mode separation is done on
		the first nCdims parameters.

	@param null_log_evidence:
		If mmodal is T, MultiNest can find multiple modes & also specify 
		which samples belong to which mode. It might be
		desirable to have separate samples & mode statistics for modes 
		with local log-evidence value greater than a
		particular value in which case nullZ should be set to that 
		value. If there isn't any particulrly interesting
		nullZ value, then nullZ should be set to a very large negative 
		number (e.g. -1.d90).
	"""

	if n_params == None:
		n_params = n_dims
	if n_clustering_params == None:
		n_clustering_params = n_params
	if wrapped_params == None:
		wrapped_params = [0] * n_params
	
	WrappedType = c_int * len(wrapped_params)
	wraps = WrappedType(*wrapped_params)
	
	prior_type = CFUNCTYPE(c_void_p, POINTER(c_double), c_int, c_int)
	loglike_type = CFUNCTYPE(c_double, POINTER(c_double), c_int, c_int)
	
	lib.run(c_int(multimodal), c_int(const_efficiency_mode), 
		c_int(n_live_points), c_double(evidence_tolerance), 
		c_double(sampling_efficiency), c_int(n_dims), c_int(n_params),
		c_int(n_clustering_params), c_int(max_modes), 
		c_int(n_iter_before_update), c_double(evidence_tolerance), 
		outputfiles_basename, c_int(seed), wraps,
		c_int(verbose), c_int(resume), 
		c_int(context), prior_type(Prior), loglike_type(LogLikelihood))

import threading

class WatchProgressThread(threading.Thread):
	"""
		Watches the progress of MultiNest.
	"""
	def __init__(self, n_params, interval_ms = 200, outputfiles_basename = "chains/1-"):
		threading.Thread.__init__(self)
		self.n_params = n_params
		self.outputfiles_basename = outputfiles_basename
		self.interval_ms = interval_ms
		"""
		This file contains the current set of live points. It has nPar+2 
		columns. The first nPar columns are the ndim
		parameter values along with the (nPar-ndim)  additional 
		parameters that are being passed by the likelihood
		routine for MultiNest to save along with the ndim parameters. 
		The nPar+1 column is the log-likelihood value &
		the last column is the node no. (used for clustering).
		"""
		self.live = "%s%s" % (self.outputfiles_basename , "phys_live.points")
		"""
		This file contains the set of rejected points. It has nPar+3 
		columns. The first nPar columns are the ndim
		parameter values along with the (nPar-ndim)  additional 
		parameters that are being passed by the likelihood
		routine for MultiNest to save along with the ndim parameters. 
		The nPar+1 column is the log-likelihood value,
		nPar+2 column is the log(prior mass) & the last column  is the 
		node no. (used for clustering).
		"""		
		self.rejected = "%s%s" % (self.outputfiles_basename , "ev.dat")
		self.running = True

	def run(self):
		import time
		while self.running:
			time.sleep(self.interval_ms / 1000.)
			try:
				print 'rejected points: ', len(file(self.rejected, 'r').readlines())
				self._plot_live()
				print 'alive points: ', len(file(self.live, 'r').readlines())
			except Exception as e:
				print e
	
	def _plot_live(self):
		import matplotlib.pyplot as plot
		import numpy
		import shutil
		x = numpy.loadtxt(self.live)
		for i in range(self.n_params):
			plot.subplot(self.n_params, 1, i)
			plot.plot(x[:,i], x[:,self.n_params], '.')
		f = "%s.pdf" % self.live
		ftmp = "%s.t.pdf" % self.live
		plot.savefig(ftmp)
		shutil.move(ftmp, f)
		print "saved to %s" % f
	


class Analyzer(object):
	def __init__(self, n_params, outputfiles_basename = "chains/1-"):
		self.outputfiles_basename = outputfiles_basename
		self.n_params = n_params
		"""
		[root].txt
		Compatable with getdist with 2+nPar columns. Columns have 
		sample probability, -2*loglikehood, samples. Sample
		probability is the sample prior mass multiplied by its 
		likelihood & normalized by the evidence.
		"""

		"""[root]post_separate.dat
		This file is only created if mmodal is set to T. Posterior 
		samples for modes with local log-evidence value
		greater than nullZ, separated by 2 blank lines. Format is the 
		same as [root].txt file."""
		self.post = "%spost_seperate.dat" % self.outputfiles_basename

		"""[root]stats.dat
		Contains the global log-evidence, its error & local log-evidence 
		with error & parameter means & standard
		deviations as well as the  best fit & MAP parameters of each of 
		the mode found with local log-evidence > nullZ.
		"""
		self.stats = "%sstats.dat" % self.outputfiles_basename

		"""
		[root]post_equal_weights.dat
		Contains the equally weighted posterior samples
		"""
		self.equal_weighted = "%spost_equal_weights.dat" % self.outputfiles_basename

	def _read_error_line(self, l):
		#print '_read_error_line', l
		name, values = l.split('  ', 1)
		name = name.strip(': ').strip()
		v, error = values.split(" +/- ")
		return name, float(v), float(error)
	def _read_error_into_dict(self, l, d):
		name, v, error = self._read_error_line(l)
		d[name.lower()] = v
		d['%s error' % name.lower()] = error
	def _read_table(self, txt, d = None, title = None):
		import numpy
		if title is None:
			title, table = txt.split("\n", 1)
		else:
			table = txt		
		header, table = txt.split("\n", 1)
		from StringIO import StringIO
		data = numpy.loadtxt(StringIO(table))
		if d is not None:
			d[title.strip().lower()] = data
		return data

	def get_stats(self):
		import re
		lines = file(self.stats).readlines()
		text = "".join(lines)
		parts = text.split("\n\n\n")
		del parts[0]
		stats = {
			'modes':[]
		}
		# Global Evidence
		self._read_error_into_dict(lines[0], stats)
		for p in parts:
			modelines = p.split("\n\n")
			mode = {
			}
			modelines1 = modelines[0].split("\n")
			# Strictly local evidence
			self._read_error_into_dict(modelines1[1], mode)
			self._read_error_into_dict(modelines1[2], mode)
			t = self._read_table(modelines[1], title = "Parameters")
			mode['mean'] = t[:,1].tolist()
			mode['sigma'] = t[:,2].tolist()
			mode['maximum'] = self._read_table(modelines[1])[:,1].tolist()
			mode['maximum a posterior'] = self._read_table(modelines[1])[:,1].tolist()
			stats['modes'].append(mode)
		return stats

if __name__ == "__main__":
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
	progress = WatchProgressThread(2)
	#progress.start()
	run(myloglike, myprior, 2)
	progress.running = False
	a = Analyzer(2)
	s = a.get_stats()

	import json
	json.dump(s, file('%s.json' % a.outputfiles_basename, 'w'), indent=2)
	print
	print "Global Evidence:\n\t%.15e +- %.15e" % ( s['global evidence'], s['global evidence error'] )
	

	"""
	create n-dimensional interpolation using scipy.interpolate.griddata
	
	create 2D marginalization plots
	add modes as 2D-ellipses; colors should add up cumulatively
	add colourmap
	add contours -- sort cubes by value; add while evidence is < Jeffrey's Factors
	
	create 1D marginalization plots
	add mode parameter estimations; as (cumulative) bars
	
	"""
	

