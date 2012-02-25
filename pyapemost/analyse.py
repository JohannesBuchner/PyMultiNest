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
	paramnames = numpy.loadtxt('params', dtype='S')[:,3]
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

