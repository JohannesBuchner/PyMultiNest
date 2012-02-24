import pymultinest
import math
import os
if not os.path.exists("chains"): os.mkdir("chains")
def show(filepath):
	""" open the output (pdf) file for the user """
	import subprocess, os
	if os.name == 'mac': subprocess.call(('open', filepath))
	elif os.name == 'nt': os.startfile(filepath)
	elif os.name == 'posix': subprocess.call(('xdg-open', filepath))

# our probability functions
# Taken from the eggbox problem.

def myprior(cube, ndim, nparams):
	#print "cube before", [cube[i] for i in range(ndim)]
	for i in range(ndim):
		cube[i] = cube[i] * 10 * math.pi
	#print "python cube after", [cube[i] for i in range(ndim)]

def myloglike(cube, ndim, nparams):
	chi = 1.
	#print "cube", [cube[i] for i in range(ndim)], cube
	for i in range(ndim):
		chi *= math.cos(cube[i] / 2.)
	#print "returning", math.pow(2. + chi, 5)
	return math.pow(2. + chi, 5)

# number of dimensions our problem has
n_params = 2

# we want to see some output while it is running
progress = pymultinest.ProgressPlotter(n_params = n_params)
progress.start()
show("chains/1-phys_live.points.pdf")
# run MultiNest
pymultinest.run(myloglike, myprior, n_params, resume = True, verbose = True, sampling_efficiency = 0.3)
# ok, done. Stop our progress watcher
progress.stop()

# lets analyse the results
a = pymultinest.Analyzer(n_params = n_params)
s = a.get_stats()

import json
json.dump(s, file('%s.json' % a.outputfiles_basename, 'w'), indent=2)
print
print "-" * 30, 'ANALYSIS', "-" * 30
print "Global Evidence:\n\t%.15e +- %.15e" % ( s['global evidence'], s['global evidence error'] )

import matplotlib.pyplot as plt
plt.clf()

# Here we will plot all the marginals and whatnot, just to show off
# You may configure the format of the output here, or in matplotlibrc
# All pymultinest does is filling in the data of the plot.

# Copy and edit this file, and play with it.

p = pymultinest.PlotMarginal(a)
for i in range(n_params):
	outfile = '%s-marginal-%d.pdf' % (a.outputfiles_basename,i)
	p.plot_marginal(i, with_ellipses = True, with_points = False, grid_points=200)
	plt.savefig(outfile, format='pdf', bbox_inches='tight')
	plt.close()
	
	outfile = '%s-mode-marginal-%d.pdf' % (a.outputfiles_basename,i)
	p.plot_modes_marginal(i, with_ellipses = True, with_points = False)
	plt.savefig(outfile, format='pdf', bbox_inches='tight')
	plt.close()
	
	outfile = '%s-mode-marginal-cumulative-%d.pdf' % (a.outputfiles_basename,i)
	p.plot_modes_marginal(i, cumulative = True, with_ellipses = True, with_points = False)
	plt.savefig(outfile, format='pdf', bbox_inches='tight')
	plt.close()
	
	for j in range(i):
		p.plot_conditional(i, j, with_ellipses = True, with_points = False)
		outfile = '%s-conditional-%d-%d.pdf' % (a.outputfiles_basename,i,j)
		plt.savefig(outfile, format='pdf', bbox_inches='tight')
		plt.close()
		show(outfile)
print "take a look at the pdf files in chains/" 



 
