"""
Module for analysing the output of APEMoST

The output files are in an easy-to-process format, so
there is no requirement to use a specific set of tools.
"""
import numpy
import scipy
import numpy
import json
import sys

def load_params():
	dtype = [('initial', 'f'), ('min', 'f'), ('max', 'f'), ('name', 'S100'), ('stepsize', 'f')]
	return numpy.loadtxt('params', dtype=dtype, ndmin=1)

def create_histogram(parameter_name, nbins=100, writeFile=True):
	"""
	Returns a histogram and some statistics about this parameter.
		
	@param writeFile: if true, write the histogram to paramname.histogram
	"""
	f = "%s-chain-0.prob.dump" % parameter_name
	values = numpy.recfromtxt(f)

	statistics = {
		'min':   float(values.min()),
		'max':   float(values.max()),
		'stdev': float(values.std()),
		'mean':  float(values.mean()),
		'median':float(numpy.median(values))
	}

	hist = scipy.histogram(values, bins=nbins, normed=True)
	histwithborders = numpy.dstack([hist[1][0:nbins], hist[1][1:nbins+1], hist[0]])
	if writeFile:
		scipy.savetxt('%s.histogram' % parameter_name, histwithborders[0], delimiter="\t")
	return histwithborders[0], statistics

def create_histograms(nbins=100, writeFile=True):
	"""
	Runs create_histogram for all parameters and returns
	a dictionary of the results
	"""
	paramnames = load_params()['name']
	return dict([(p, create_histogram(p, nbins=nbins, writeFile=writeFile)) for p in paramnames])

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
	
	print "Model probability ln(p(D|M, I)): [about 10^%.0f] %.5f" % (logprob / scipy.log(10), logprob)
	print ("""
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
	""" % tuple(['s']*10)) % limits

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
		value = numpy.loadtxt(f).mean() / beta[i]
		logprob = logprob + value * (beta[i] - previous_beta)
		previous_beta = beta[i]
	
	if show is True:
		print_model_probability(logprob)
	return logprob

import numpy
import matplotlib
import matplotlib.pyplot as plt

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
	
	def plot(self):
		# load data
		values = []
		print "loading chains ..."
		for p in self.params:
			f = "%s-chain-0.prob.dump" % p['name']
			print "	loading chain %s" % f
			v = numpy.genfromtxt(f, skip_footer=1, dtype='f')
			values.append(v)
		print "loading chains finished."
		nvalues = min(map(len, values))
		values = map(lambda v: v[-self.nlast:nvalues], values)
		
		for p1,v1,i in zip(self.params, values, range(len(self.params))):
			self.marginal_plot(p1, v1)
			for p2,v2 in zip(self.params[:i], values[:i]):
				self.conditional_plot(p1, v1, p2, v2)
	
class VisitedPlotter(VisitedAnalyser):
	"""
	Produces individual plots of visited points for each parameter,
	and each parameter pair.
	
	@param outputfiles_basename: prefix of output files
	
	The output files are named chain0-paramname1-paramname2.pdf and chain0-paramname.pdf
	"""
	def __init__(self, outputfiles_basename = "", nlast = 0):
		VisitedAnalyser.__init__(self, nlast = 0)
		self.outputfiles_basename = outputfiles_basename
	
	def conditional_plot(self, param1, values1, param2, values2):
		self.conditional_plot_before(param1, values1, param2, values2)
		plt.plot(values1, values2, '+', color='black', markersize=2, alpha=0.5)
		self.conditional_plot_after(param1, values1, param2, values2)
	
	def conditional_plot_before(self, param1, values1, param2, values2):
		names = (param1['name'],param2['name'])
		print "creating conditional plot of %s vs %s" % names
		plt.figure(figsize=(5,5))
		plt.title("%s vs %s" % names)
	def conditional_plot_after(self, param1, values1, param2, values2):
		names = (param1['name'],param2['name'])
		plt.savefig(self.outputfiles_basename + "chain0-%s-%s.pdf" % names)
	
	def marginal_plot(self, param, values):
		self.marginal_plot_before(param, values)
		plt.plot(values, color='gray')
		self.marginal_plot_after(param, values)
	def marginal_plot_before(self, param, values):
		name = param['name']
		print "creating marginal plot of %s" % name
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
	def plot(self):
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
		
		plt.figure(figsize=(5*self.nparams,5*self.nparams))
		VisitedPlotter.plot(self)
		print "saving output ..."
		plt.savefig(self.outputfiles_basename + "chain0.pdf")
		print "saving output done"
	
	def choose_plot(self, i, j):
		plt.subplot(self.nparams, self.nparams, self.nparams * j + i + 1)
	
	def conditional_plot_before(self, param1, values1, param2, values2):
		names = [param1['name'],param2['name']]
		i, j = map(self.paramnames.index, names)
		self.choose_plot(i, j)
		names = (param1['name'],param2['name'])
		print "creating conditional plot of %s vs %s" % names
		plt.xlabel(param1['name'])
		plt.ylabel(param2['name'])
	def conditional_plot_after(self, param1, values1, param2, values2):
		pass
	
	def marginal_plot_before(self, param, values):
		name = param['name']
		i = self.paramnames.index(name)
		thisplot = self.choose_plot(i, i)
		print "creating marginal plot of %s" % name
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
		self.update_plot(thisplot, "progress", xrange(len(values)), values)
		plt.xlabel("iteration")
		plt.ylabel(name)
		
	def plot(self):
		VisitedAnalyser.plot(self)
		plt.show()
