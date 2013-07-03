from __future__ import absolute_import, unicode_literals, print_function
import scipy
import numpy
import math
import sys
import matplotlib.pyplot as plt
import pyapemost

def show(filepath):
	""" open the output (pdf) file for the user """
	import subprocess, os
	if os.name == 'mac': subprocess.call(('open', filepath))
	elif os.name == 'nt': os.startfile(filepath)
	elif os.name == 'posix': subprocess.call(('xdg-open', filepath))

# we will try to recover these parameters
A     = 6.13
omega = 1
phase = 0.314
sigma = 10
offset = 0

# generate some data
x = numpy.linspace(0, 10, 1000)
y = scipy.random.normal(A * numpy.sin((x * omega + phase) * 2 * math.pi) + offset, sigma)
numpy.savetxt("data", numpy.vstack((x,y)).transpose())
plt.title("Problem: Finding a sine signal in noisy data")
plt.suptitle("amplitude: %f, phase: %f, frequency: %f, offset: %f" % (A, phase, omega, offset))
plt.plot(x, A * numpy.sin((x * omega + phase) * 2 * math.pi), '--r', label="true sin")
plt.plot(x, y, 'ob', label="observed data")
plt.legend(loc='best')
plt.savefig("data.pdf")
plt.clf()
show("data.pdf")

f = file("params",'w')
f.write("0	0	10	amplitude	-1\n")
f.write("0	0	40	frequency	-1\n")
f.write("0	0	1	phase    	-1\n")
f.write("0	-10	10	offset   	-1\n")
f.close()

# our probability functions
# A simple chi-square of a sine model.

def myprior(m, oldparams):
	return 0.

def myloglike(mp, oldparams):
	m = mp.contents
	amplitude, frequency, phase, offset = [m.params.contents.data[i] for i in range(m.n_params)]
	
	y_model = amplitude * numpy.sin(2.0 * numpy.pi * (frequency * x + phase)) + offset
	squares = (y_model - y)**2
	prob = squares.sum() / (-2. * sigma**2)
	#print amplitude, frequency, phase, offset, prob
	return prob

pyapemost.set_function(myloglike, myprior)

if len(sys.argv) < 2:
	#print "SYNAPSIS: %s [calibrate|run|analyse]" % sys.argv[0]
	cmd = ["calibrate", "run", "analyse"]
else:
	cmd = sys.argv[1:]

if "calibrate" in cmd:
	pyapemost.calibrate()
# run APEMoST
if "run" in cmd:
	w = pyapemost.watch.ProgressWatcher()
	w.start()
	pyapemost.run(max_iterations = 10000, append = False)
	w.stop()

# lets analyse the results
if "analyse" in cmd:
	plotter = pyapemost.analyse.VisitedAllPlotter()
	plotter.plot()
	show("chain0.pdf")
	histograms = pyapemost.create_histograms()
	i = 1
	plt.clf()
	plt.figure(figsize=(7, 4 * len(histograms)))
	for k,(v,stats) in histograms.items():
		plt.subplot(len(histograms), 1, i)
		plt.plot(v[:,0], v[:,2], ls='steps--', label=k)
		plt.legend()
		print(k, stats)
		i = i + 1
	plt.savefig("marginals.pdf")
	show("marginals.pdf")
	
	print(pyapemost.model_probability())


