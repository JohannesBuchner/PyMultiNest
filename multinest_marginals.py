#!/usr/bin/env python
__doc__ = """
Script that does default visualizations (marginal plots, 1-d and 2-d).

Author: Johannes Buchner (C) 2013
"""

import numpy
from numpy import exp, log
import matplotlib.pyplot as plt
import sys, os
import json
import pymultinest

if len(sys.argv) != 2:
	sys.stderr.write("""SYNOPSIS: %s <output-root> 

	output-root: 	Where the output of a MultiNest run has been written to. 
	            	Example: chains/1-
%s""" % (sys.argv[0], __doc__))
	sys.exit(1)

prefix = sys.argv[1]

parameters = json.load(file(prefix + 'params.json'))
n_params = len(parameters)

a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
s = a.get_stats()

p = pymultinest.PlotMarginal(a)

values = a.get_equal_weighted_posterior()
assert n_params == len(s['marginals'])
modes = s['modes']

dim2 = os.environ.get('D', '1' if n_params > 20 else '2') == '2'
nbins = 100 if n_params < 3 else 20
if dim2:
	plt.figure(figsize=(5*n_params, 5*n_params))
	for i in range(n_params):
		plt.subplot(n_params, n_params, i + 1)
		plt.xlabel(parameters[i])
	
		m = s['marginals'][i]
		plt.xlim(m['5sigma'])
	
		oldax = plt.gca()
		x,w,patches = oldax.hist(values[:,i], bins=nbins, edgecolor='grey', color='grey', histtype='stepfilled', alpha=0.2)
		oldax.set_ylim(0, x.max())
	
		newax = plt.gcf().add_axes(oldax.get_position(), sharex=oldax, frameon=False)
		p.plot_marginal(i, ls='-', color='blue', linewidth=3)
		newax.set_ylim(0, 1)
	
		ylim = newax.get_ylim()
		y = ylim[0] + 0.05*(ylim[1] - ylim[0])
		center = m['median']
		low1, high1 = m['1sigma']
		print center, low1, high1
		newax.errorbar(x=center, y=y,
			xerr=numpy.transpose([[center - low1, high1 - center]]), 
			color='blue', linewidth=2, marker='s')
		oldax.set_yticks([])
		#newax.set_yticks([])
		newax.set_ylabel("Probability")
		ylim = oldax.get_ylim()
		newax.set_xlim(m['5sigma'])
		oldax.set_xlim(m['5sigma'])
		#plt.close()
	
		for j in range(i):
			plt.subplot(n_params, n_params, n_params * (j + 1) + i + 1)
			p.plot_conditional(i, j, bins=20, cmap = plt.cm.gray_r)
			for m in modes:
				plt.errorbar(x=m['mean'][i], y=m['mean'][j], xerr=m['sigma'][i], yerr=m['sigma'][j])
			plt.xlabel(parameters[i])
			plt.ylabel(parameters[j])
			#plt.savefig('cond_%s_%s.pdf' % (params[i], params[j]), bbox_tight=True)
			#plt.close()

	plt.savefig(prefix + 'marg.pdf')
	plt.savefig(prefix + 'marg.png')
	plt.close()
else:
	from matplotlib.backends.backend_pdf import PdfPages
	print '1dimensional only. Set the D environment variable D=2 to force'
	print '2d marginal plots.'
	pp = PdfPages(prefix + 'marg1d.pdf')
	
	for i in range(n_params):
		plt.figure(figsize=(5, 5))
		plt.xlabel(parameters[i])
		
		m = s['marginals'][i]
		plt.xlim(m['5sigma'])
	
		oldax = plt.gca()
		x,w,patches = oldax.hist(values[:,i], bins=20, edgecolor='grey', color='grey', histtype='stepfilled', alpha=0.2)
		oldax.set_ylim(0, x.max())
	
		newax = plt.gcf().add_axes(oldax.get_position(), sharex=oldax, frameon=False)
		p.plot_marginal(i, ls='-', color='blue', linewidth=3)
		newax.set_ylim(0, 1)
	
		ylim = newax.get_ylim()
		y = ylim[0] + 0.05*(ylim[1] - ylim[0])
		center = m['median']
		low1, high1 = m['1sigma']
		print center, low1, high1
		newax.errorbar(x=center, y=y,
			xerr=numpy.transpose([[center - low1, high1 - center]]), 
			color='blue', linewidth=2, marker='s')
		oldax.set_yticks([])
		newax.set_ylabel("Probability")
		ylim = oldax.get_ylim()
		newax.set_xlim(m['5sigma'])
		oldax.set_xlim(m['5sigma'])
		plt.savefig(pp, format='pdf', bbox_inches='tight')
		plt.close()
	pp.close()



