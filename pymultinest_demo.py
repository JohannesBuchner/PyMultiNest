from __future__ import absolute_import, unicode_literals, print_function
import pymultinest as pymn
import math
import os, threading, subprocess
if not os.path.exists("chains"): os.mkdir("chains")
def show(filepath): 
	""" open the output (pdf) file for the user """
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
parameters = ["x", "y"]
n_params = len(parameters)

# we want to see some output while it is running
progress = pymn.ProgressPlotter(n_params = n_params); progress.start()
threading.Timer(2, show, ["chains/1-phys_live.points.pdf"]).start() # delayed opening
# run MultiNest
pymn.run(myloglike, myprior, n_params, importance_nested_sampling = False, resume = True, verbose = True, sampling_efficiency = 0.3)
# ok, done. Stop our progress watcher
progress.stop()

# lets analyse the results
a = pymn.Analyzer(n_params = n_params)
s = a.get_stats()

import json
#json.dump(s, file('%s.json' % a.outputfiles_basename, 'w'), indent=2)
with open('%s.json' % a.outputfiles_basename, mode='w') as file:
	json.dump(s, file, indent=2)
	print()
	print("-" * 30, 'ANALYSIS', "-" * 30)
	print("Global Evidence:\n\t%.15e +- %.15e" % ( s['nested sampling global log-evidence'], s['nested sampling global log-evidence error'] ))

import matplotlib.pyplot as plt
plt.clf()

# Here we will plot all the marginals and whatnot, just to show off
# You may configure the format of the output here, or in matplotlibrc
# All pymultinest does is filling in the data of the plot.

# Copy and edit this file, and play with it.

p = pymn.PlotMarginalModes(a)
plt.figure(figsize=(5*n_params, 5*n_params))
#plt.subplots_adjust(wspace=0, hspace=0)
for i in range(n_params):
	plt.subplot(n_params, n_params, n_params * i + i + 1)
	p.plot_marginal(i, with_ellipses = True, with_points = False, grid_points=50)
	plt.ylabel("Probability")
	plt.xlabel(parameters[i])
	
	for j in range(i):
		plt.subplot(n_params, n_params, n_params * j + i + 1)
		#plt.subplots_adjust(left=0, bottom=0, right=0, top=0, wspace=0, hspace=0)
		p.plot_conditional(i, j, with_ellipses = False, with_points = True, grid_points=30)
		plt.xlabel(parameters[i])
		plt.ylabel(parameters[j])

plt.savefig("chains/marginals_multinest.pdf") #, bbox_inches='tight')
show("chains/marginals_multinest.pdf")
for i in range(n_params):
	outfile = '%s-mode-marginal-%d.pdf' % (a.outputfiles_basename,i)
	p.plot_modes_marginal(i, with_ellipses = True, with_points = False)
	plt.ylabel("Probability")
	plt.xlabel(parameters[i])
	plt.savefig(outfile, format='pdf', bbox_inches='tight')
	plt.close()
	
	outfile = '%s-mode-marginal-cumulative-%d.pdf' % (a.outputfiles_basename,i)
	p.plot_modes_marginal(i, cumulative = True, with_ellipses = True, with_points = False)
	plt.ylabel("Cumulative probability")
	plt.xlabel(parameters[i])
	plt.savefig(outfile, format='pdf', bbox_inches='tight')
	plt.close()

print("Take a look at the pdf files in chains/") 



 
