from __future__ import print_function
import numpy
from numpy import cos, pi


import pypolychord

def loglikelihood(params):
	return (numpy.prod(cos(params / 2)) + 2)**5

def prior(pars):
	return (4 + 1) * pi * pars

if __name__ == "__main__":
	ndim = 2
	import os
	if not os.path.exists('chains'): os.mkdir('chains')
	if not os.path.exists('chains/clusters'): os.mkdir('chains/clusters')
	pypolychord.run(loglikelihood, prior, ndim, n_live=500, n_chords=1, output_basename='chains/eggbox-')


