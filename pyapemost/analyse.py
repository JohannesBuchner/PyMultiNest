"""
Module for analysing the output of APEMoST

The output files are in an easy-to-process format, so
there is no requirement to use a specific set of tools.
"""
from __future__ import absolute_import, unicode_literals, print_function
import numpy
import scipy, scipy.stats
import numpy
import json
import sys, os

nevery = 1

def load_params():
	dtype = [('initial', 'f'), ('min', 'f'), ('max', 'f'), ('name', 'S100'), ('stepsize', 'f')]
	return numpy.loadtxt('params', dtype=dtype, ndmin=1)

import scipy.ndimage.filters
def smoothen_histogram(bins):
	bins[:,2] = scipy.ndimage.filters.median_filter(bins[:,2], size=10)
	return bins

def truncate_histogram(bins, sigma):
	#print 'bins', bins[:,2]
	bins_sorted = numpy.asarray(sorted(bins[:,2].reshape((-1,)), reverse=True))
	#print 'bins_sorted', bins_sorted
	bins_cum = bins_sorted.cumsum()
	bins_tot = bins_sorted.sum()
	#print 'bins_cum', bins_cum
	
	limit = scipy.stats.norm.cdf(sigma) - scipy.stats.norm.cdf(-sigma)
	badbins = bins_cum > limit * bins_tot
	
	# find next-best
	#print "limit for %f is %f: %d/%d bins" % (sigma, limit, badbins.sum(), len(bins))
	if badbins.sum() != 0:
		minvalue = bins_sorted[badbins][0]
		
		#print 'good bins', bins[bins[:,2] >= minvalue]
		bins[bins[:,2] < minvalue,2] = 0.
		return bins # [bins[:,2] >= minvalue]
	else:
		return bins

def create_histogram(parameter_name, nbins=100, writeFile=True, skipfirst=0, truncate=False, smooth=False):
	"""
	Returns a histogram and some statistics about this parameter.
		
	@param writeFile: if true, write the histogram to paramname.histogram
	"""
	f = "%s-chain-0.prob.dump" % parameter_name
	values = numpy.recfromtxt(f)[skipfirst::nevery]

	statistics = {
		'min':   float(values.min()),
		'max':   float(values.max()),
		'stdev': float(values.std()),
		'mean':  float(values.mean()),
		'median':float(numpy.median(values)),
		'q1':    float(scipy.stats.scoreatpercentile(values, 25)),
		'q3':    float(scipy.stats.scoreatpercentile(values, 75)),
		'p5':    float(scipy.stats.scoreatpercentile(values, 5)),
		'p95':    float(scipy.stats.scoreatpercentile(values, 95)),
	}
	
	hist = scipy.histogram(values, bins = nbins if not smooth else nbins*10, normed=True)
	histwithborders = numpy.dstack([hist[1][0:nbins], hist[1][1:nbins+1], hist[0]])
	if writeFile:
		scipy.savetxt('%s.histogram' % parameter_name, histwithborders[0], delimiter="\t")
	return histwithborders[0], statistics


def create_histograms(**kwargs):
	"""
	Runs create_histogram for all parameters and returns
	a dictionary of the results
	"""
	paramnames = load_params()['name']
	return dict([(p, create_histogram(p, **kwargs)) for p in paramnames])

def print_model_probability(logprob):
	"""
	Gives a nice overview of the model probability, allowing
	the practitioner to compare this model's probability to others
	"""
	prob = scipy.exp(logprob)

	limits = {
		'eq'         :logprob, 
		'barely'     :logprob - scipy.log(3),
		'substantial':logprob - scipy.log(10),
		'strong'     :logprob - scipy.log(30),
		'very strong':logprob - scipy.log(100)
	}
	for i in limits:
		limits[i] = "%5.1f" % limits[i]
		limits[i] = " " * (7 - len(limits[i])) + limits[i]
	
	print("Model probability ln(p(D|M, I)): [about 10^%.0f] %.5f" % (logprob / scipy.log(10), logprob))
	print(("""
	Table to compare support against other models (Jeffrey):

	   other model     |
	   ln(p(D|M,I))    | supporting evidence for this model
	 ------------------+-------------------------------------
		  >%%(eq)%s   Negative (supports other model)
	 %%(eq)%s ..%%(barely)%s   Barely worth mentioning
	 %%(barely)%s ..%%(substantial)%s   Substantial
	 %%(substantial)%s ..%%(strong)%s   Strong
	 %%(strong)%s ..%%(very strong)%s   Very strong
		  <%%(very strong)%s   Decisive

	be careful.
	""" % tuple(['s']*10)) % limits)

def model_probability(show=True):
	"""
	Calculate the model probability.
	
	@param show: If true, print_model_probability is called
	"""
	calibration = numpy.loadtxt('calibration_results', ndmin=2)
	nchains = len(calibration)
	beta = calibration[:,0]

	previous_beta = 0.
	logprob = 0.
	for i in range(nchains - 1, -1, -1):
		f = 'prob-chain%d.dump' % i
		# this file has 2 columns: with prior and *without*.
		value = numpy.loadtxt(f)[:,1].mean() / beta[i]
		logprob = logprob + value * (beta[i] - previous_beta)
		previous_beta = beta[i]
	
	if show is True:
		print_model_probability(logprob)
	return logprob

import numpy
import matplotlib
import matplotlib.pyplot as plt

verbose = 1

class VisitedAnalyser(object):
	"""
	Internal, abstract class that has the visited points of chain 0 loaded. 
	
	You need to overwrite::
	
		marginal_plot(param1, values1)
		conditional_plot(param1, values1, param2, values2)
		
	@param nlast: if 0, use all points; otherwise only the last *nlast* ones.
	"""
	def __init__(self, nlast = 0):
		self.params = load_params()
		self.nlast = nlast
	
	def load(self):
		# load data
		values = []
		if verbose: print()
		if verbose: print("visualization: loading chains ...")
		f = "prob-chain0.dump"
		if not os.path.exists(f):
			raise Exception("visualization: chains not available yet.")
		try:
			# I think the first column is the probabilities, the second is without prior
			probabilities = numpy.genfromtxt(f, skip_footer=1, dtype='f')[:,0]
		except Exception as e:
			raise Exception("visualization: chains couldn't be loaded; perhaps no data yet: " + str(e))
		for p in self.params:
			f = "%s-chain-0.prob.dump" % p['name']
			if verbose: print("	loading chain %s" % f)
			if not os.path.exists(f):
				raise Exception("visualization: chains not available yet.")
			try:
				v = numpy.genfromtxt(f, skip_footer=1, dtype='f')
			except Exception as e:
				raise Exception("visualization: chains couldn't be loaded; perhaps no data yet: " + str(e))
			values.append(v)
		nvalues = min(map(len, values))
		if verbose: print("visualization: loading chains finished; %d values" % nvalues)
		self.values = [v[:nvalues][-self.nlast::nevery] for v in values]
		self.probabilities = probabilities[:nvalues][-self.nlast::nevery]
	
	def plot_only(self):
		for p1,v1,i in zip(self.params, self.values, range(len(self.params))):
			self.marginal_plot(p1, v1, self.probabilities)
			for p2,v2 in zip(self.params[:i], self.values[:i]):
				self.conditional_plot(p1, v1, p2, v2, self.probabilities)
	def after_plot(self):
		pass
	def before_plot(self):
		pass
	def plot(self):
		self.load()
		self.before_plot()
		self.plot_only()
		self.after_plot()
		del self.values
		del self.probabilities
	
class VisitedPlotter(VisitedAnalyser):
	"""
	Produces individual plots of visited points for each parameter,
	and each parameter pair.
	
	@param outputfiles_basename: prefix of output files
	
	The output files are named chain0-paramname1-paramname2.pdf and chain0-paramname.pdf
	"""
	def __init__(self, outputfiles_basename = "", nlast = 0):
		VisitedAnalyser.__init__(self, nlast = nlast)
		self.outputfiles_basename = outputfiles_basename
	
	def conditional_plot(self, param1, values1, param2, values2, probabilities, points = True, contours = False, best = True, good = True, sigmas = [1.,3.], contourargs = {}):
		self.conditional_plot_before(param1, values1, param2, values2)
		if best:
			best = probabilities.argmax()
		if good:
			good = probabilities > (probabilities[best] - 0.5)
		else:
			good = numpy.zeros_like(probabilities)
		
		CS = None
		sigmas = numpy.asarray(sigmas)
		
		if points:
			plt.plot(values1[-good], values2[-good], '+', color='#444444', markersize=2, alpha=0.5)
			plt.plot(values1[good], values2[good], 'x', color='#22FF22', markersize=2, alpha=0.5, label='good')
			if best:
				plt.plot([values1[best]], [values2[best]], 'o', color='red', markersize=2, alpha=0.5, label='best')
		
		if contours:
			H, xedges, yedges = numpy.histogram2d(values1, values2, bins=int(len(values1)**0.5 / 4), normed=True)
			#print 'edges', xedges, yedges
			
			# using the binning results as the value of the center of the histogram
			# would be an inaccurate visualization of the borders
			# instead, we make 4 points (one for each corner)
			
			#X1, Y1 = numpy.meshgrid(
				#(xedges[:-1]*9 + xedges[1:])/10., 
				#(yedges[:-1]*9 + yedges[1:])/10.)
			Z1 = H / H.max()
			
			X = []
			Y = []
			Z = []
			for i in range(len(xedges) - 1):
				xlow = xedges[i]
				xhigh = xedges[i+1]
				X_row_low  = [(xlow*9 + xhigh)/10.]*((len(yedges)-1)*2)
				X_row_high = [(xlow + 9*xhigh)/10.]*((len(yedges)-1)*2)
				Y_row = []
				Z_row = []
				
				for j in range(len(yedges)-1):
					ylow = yedges[j]
					yhigh = yedges[j+1]
					#print ylow, yhigh
					Y_row.append((ylow*9 + yhigh)/10.)
					Y_row.append((ylow + 9*yhigh)/10.)
					z = Z1[i,j]
					Z_row += [z]*2
				#print len(X_row_low), len(X_row_high), len(yedges)
				X.append(X_row_low)
				X.append(X_row_high)
				Y.append(Y_row)
				Y.append(Y_row)
				Z.append(Z_row)
				Z.append(Z_row)
			
			X = numpy.array(X)
			Y = numpy.array(Y)
			Z = numpy.array(Z)
			#print X.shape, Y.shape, Z.shape
			
			ndim = len(self.params)
			
			# order bins by value
			#print 'Z', Z
			z_sorted = numpy.asarray(sorted(Z1.reshape((-1,)), reverse=True))
			#print 'z_sorted', z_sorted
			# adding up largest bins first
			z_cum = z_sorted.cumsum()
			z_tot = z_sorted.sum()
			#print 'z_cum', z_cum
			
			z_mins = []
			contour_sigmas = []
			colors = []
			for s in sigmas:
				limit = scipy.stats.norm.cdf(s) - scipy.stats.norm.cdf(-s)
				badbins = z_cum > limit * z_tot
				# find next-best
				#print "limit for %f is %f: %d bins" % (s, limit, badbins.sum())
				if badbins.sum() != 0:
					z_min = z_sorted[badbins][0]
					z_mins.append(z_min)
					contour_sigmas.append('%.0f sigma' % s)
					if s == 1: colors.append("#dddddd")
					if s == 3: colors.append("#cccccc")
			
			#print 'z_mins', z_mins
			CS = plt.contour(X, Y, Z, z_mins, labels=contour_sigmas, **contourargs)
			plt.clabel(CS, fmt=dict(zip(z_mins, contour_sigmas)), inline=1, fontsize=10)
			if len(z_mins) > 1:
				#print [z_mins[-1], z_mins[0]]
				plt.contourf(X, Y, Z, z_mins[::-1] + [Z1.max()+1], colors=colors, alpha=0.3)
			
			## 1, 3 and 5 sigma equivalents by FWHM
			##lines = 2 * (2 * ndim * numpy.log([1,3,5])**0.5
			#lines = numpy.exp(-1/2. * sigmas)
			#CS = plt.contour(X, Y, Z, lines, **contourargs)
		
		self.conditional_plot_after(param1, values1, param2, values2)
		return CS
	
	def conditional_plot_before(self, param1, values1, param2, values2):
		names = (param1['name'],param2['name'])
		if verbose: print("visualization: creating conditional plot of %s vs %s" % names)
		plt.figure(figsize=(5,5))
		plt.title("%s vs %s" % names)
	def conditional_plot_after(self, param1, values1, param2, values2):
		names = (param1['name'],param2['name'])
		plt.savefig(self.outputfiles_basename + "chain0-%s-%s.pdf" % names)
	
	def marginal_plot(self, param, values, probabilities):
		self.marginal_plot_before(param, values)
		best = probabilities.argmax()
		good = probabilities > (probabilities[best] - 0.5)
		i = numpy.arange(len(probabilities))
		plt.plot(i[-good], values[-good], '+', color='gray')
		plt.plot(i[good], values[good], 'x', color='#22FF22')
		plt.plot([i[best]], [values[best]], 'o', color='red')
		
		self.marginal_plot_after(param, values)
	def marginal_plot_before(self, param, values):
		name = param['name']
		if verbose: print("visualization: creating marginal plot of %s" % name)
		plt.figure(figsize=(5,5))
		plt.title("%s" % name)
	def marginal_plot_after(self, param, values):
		name = param['name']
		plt.savefig(self.outputfiles_basename + "chain0-%s.pdf" % name)
		plt.clf()

class VisitedAllPlotter(VisitedPlotter):
	"""
	Extends VisitedPlotter to put all of those plots into one file (named chain0.pdf).
	
	@see VisitedPlotter
	"""
	def before_plot(self):
		self.paramnames = list(self.params['name'])
		self.nparams = len(self.params)
		self.plots = {}
		for i, p1 in zip(range(self.nparams), self.params):
			self.choose_plot(i, i)
			plt.title("%s" % p1['name'])
			
			for j, p2 in zip(range(i), self.params):
				self.choose_plot(i, j)
				names = (str(p1['name']), str(p2['name']))
				plt.title("%s vs %s" % names)
			for j, p2 in zip(range(i), self.params):
				self.choose_plot(j, i)
				names = (str(p1['name']), str(p2['name']))
				plt.title("%s vs %s" % names)
		self.before_plot_newfigure()
	def before_plot_newfigure(self):
		plt.figure(figsize=(5*self.nparams,5*self.nparams))
	
	def after_plot(self):
		if verbose: print("visualization: saving output ...")
		plt.savefig(self.outputfiles_basename + "chain0.pdf")
		plt.close()
		if verbose: print("visualization: saving output done")
	
	def choose_plot(self, i, j):
		plt.subplot(self.nparams, self.nparams, self.nparams * j + i + 1)
	
	def conditional_plot_before(self, param1, values1, param2, values2):
		names = [param1['name'],param2['name']]
		i, j = map(self.paramnames.index, names)
		self.choose_plot(i, j)
		names = (param1['name'],param2['name'])
		if verbose: print("visualization: creating conditional plot of %s vs %s" % names)
		plt.xlabel(param1['name'])
		plt.ylabel(param2['name'])
	def conditional_plot(self, param1, values1, param2, values2, probabilities, **kwargs):
		VisitedPlotter.conditional_plot(self, param1, values1, param2, values2, points = True, contours = False, probabilities = probabilities, **kwargs)
		VisitedPlotter.conditional_plot(self, param2, values2, param1, values1, points = False, contours = True, probabilities = probabilities, **kwargs)
	
	def conditional_plot_after(self, param1, values1, param2, values2):
		pass
	
	def marginal_plot_before(self, param, values):
		name = param['name']
		i = self.paramnames.index(name)
		thisplot = self.choose_plot(i, i)
		if verbose: print("visualization: creating marginal plot of %s" % name)
		plt.xlabel("iteration")
		plt.ylabel(name)
	def marginal_plot_after(self, param, values):
		pass
	

class VisitedWindow(VisitedAnalyser):
	"""
	Not implemented / tested yet.
	
	Should give a window to watch the progress of the Markov Chain,
	similar to **VisitedAllPlotter** but continuously updating.
	"""
	def __init__(self):
		plt.clf()
		plt.ion()
		VisitedAnalyser.__init__(self)
		self.paramnames = list(self.params['name'])
		self.nparams = len(self.params)
		self.plots = {}
		plt.figure(figsize=(3*self.nparams,3*self.nparams))
		for i, p1 in zip(range(self.nparams), self.params):
			self.plots[i] = {i:{}}
			self.choose_plot(i, i)
			plt.title("%s" % p1['name'])
			
			for j, p2 in zip(range(i), self.params):
				self.plots[i][j] = {}
				self.choose_plot(i, j)
				names = (str(p1['name']), str(p2['name']))
				plt.title("%s vs %s" % names)
	
	def choose_plot(self, i, j):
		plt.subplot(self.nparams, self.nparams, self.nparams * j + i + 1)
		return self.plots[i][j]
	def update_plot(self, plot, name, x, y):
		if name not in plot or True:
			plot[name], = plt.plot(x, y, '+', markersize=2, alpha=0.5, label=name)
		else:
			plot[name].set_xdata(x)
			plot[name].set_ydata(y)
		plt.draw()
		plt.show()
	
	def conditional_plot(self, param1, values1, param2, values2):
		names = [param1['name'],param2['name']]
		i, j = map(self.paramnames.index, names)
		thisplot = self.choose_plot(i, j)
		
		self.update_plot(thisplot, "newest", values1[-50:], values2[-50:])
		plt.xlabel(param1['name'])
		plt.ylabel(param2['name'])
		
	def marginal_plot(self, param, values):
		name = param['name']
		i = self.paramnames.index(name)
		thisplot = self.choose_plot(i, i)
		self.update_plot(thisplot, "progress", range(len(values)), values)
		plt.xlabel("iteration")
		plt.ylabel(name)
		
	def after_plot(self):
		plt.show()
