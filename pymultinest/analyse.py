import numpy
from StringIO import StringIO
import re

class Analyzer(object):
	"""
		Class for accessing the output of MultiNest.
		
		After the run, this class provides a mean of reading the result
		in a structured way.
		
	"""
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
		self.data_file = "%s.txt" % self.outputfiles_basename
		print '  analysing data from %s' % self.data_file

		"""[root]post_separate.dat
		This file is only created if mmodal is set to T. Posterior 
		samples for modes with local log-evidence value
		greater than nullZ, separated by 2 blank lines. Format is the 
		same as [root].txt file."""
		self.post_file = "%spost_seperate.dat" % self.outputfiles_basename

		"""[root]stats.dat
		Contains the global log-evidence, its error & local log-evidence 
		with error & parameter means & standard
		deviations as well as the  best fit & MAP parameters of each of 
		the mode found with local log-evidence > nullZ.
		"""
		self.stats_file = "%sstats.dat" % self.outputfiles_basename

		"""
		[root]post_equal_weights.dat
		Contains the equally weighted posterior samples
		"""
		self.equal_weighted_file = "%spost_equal_weights.dat" % self.outputfiles_basename
	
	def get_data(self):
		"""
			fetches self.data_file
		"""
		if not hasattr(self, 'data'):
			self.data = numpy.loadtxt(self.data_file, ndmin=2)
		return self.data
	def get_equal_weighted_posterior(self):
		"""
			fetches self.data_file
		"""
		if not hasattr(self, 'equal_weighted_posterior'):
			self.equal_weighted_posterior = numpy.loadtxt(self.equal_weighted_file, ndmin=2)
		return self.equal_weighted_posterior
	
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
		if title is None:
			title, table = txt.split("\n", 1)
		else:
			table = txt		
		header, table = txt.split("\n", 1)
		data = numpy.loadtxt(StringIO(table), ndmin=2)
		if d is not None:
			d[title.strip().lower()] = data
		return data
	
	def get_stats(self):
		posterior = self.get_data()
		
		stats = []
		
		for i in range(2, posterior.shape[1]):
			b = zip(posterior[:,0], posterior[:,i])
			b.sort(key=lambda x: x[1])
			b = numpy.array(b)
			b[:,0] = b[:,0].cumsum()
			sig5 = 0.5 + 0.9999994 / 2.
			sig3 = 0.5 + 0.9973 / 2.
			sig2 = 0.5 + 0.95 / 2.
			sig1 = 0.5 + 0.6826 / 2.
			low1  = b[b[:,0] < 1 - sig1][-1][1]
			high1 = b[b[:,0] > sig1][0][1]
			low2  = b[b[:,0] < 1 - sig2][-1][1]
			high2 = b[b[:,0] > sig2][0][1]
			try:
				low3  = b[b[:,0] <= 1 - sig3][-1][1]
			except:
				low3 = b[0][1]
			try:
				high3 = b[b[:,0] >= sig3][0][1]
			except:
				high3 = b[-1][1]
			try:
				low5  = b[b[:,0] <= 1 - sig5][-1][1]
			except:
				low5 = b[0][1]
			try:
				high5 = b[b[:,0] >= sig5][0][1]
			except:
				high5 = b[-1][1]
			median = b[b[:,0] >= 0.5][0][1]
			q1 = b[b[:,0] >= 0.75][0][1]
			q3 = b[b[:,0] <= 0.25][-1][1]
			q99 = b[b[:,0] >= 0.99][0][1]
			q01 = b[b[:,0] <= 0.01][-1][1]
			q90 = b[b[:,0] >= 0.9][0][1]
			q10 = b[b[:,0] <= 0.1][-1][1]
			
			stats.append({
				'median': median,
				'sigma': (high1 - low1) / 2.,
				'1sigma': [low1, high1],
				'2sigma': [low2, high2],
				'3sigma': [low3, high3],
				'5sigma': [low5, high5],
				'q75%': q1,
				'q25%': q3,
				'q99%': q99,
				'q01%': q01,
				'q90%': q99,
				'q10%': q10,
			})
		
		mode_stats = self.get_mode_stats()
		mode_stats['marginals'] = stats 
		return mode_stats
	
	def get_mode_stats(self):
		"""
			information about the modes found:
			mean, sigma, maximum a posterior in each dimension
		"""
		lines = file(self.stats_file).readlines()
		text = "".join(lines)
		parts = text.split("\n\n\n")
		del parts[0]
		stats = {'modes':[]}
		
		# Global Evidence
		self._read_error_into_dict(lines[0], stats)
		i = 0
		for p in parts:
			modelines = p.split("\n\n")
			mode = {
				'index':i
			}
			i = i + 1
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

	def get_best_fit(self):
		lastrow = self.get_data()[-1]
		return {'log_likelihood': float(-0.5 * lastrow[1]), 
			'parameters': list(lastrow[2:])}


