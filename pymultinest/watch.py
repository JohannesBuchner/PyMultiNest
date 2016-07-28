from __future__ import absolute_import, unicode_literals, print_function
import threading

class ProgressWatcher(threading.Thread):
	"""
		Abstract class for watching the progress of MultiNest.
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
	
	def stop(self):
		self.running = False

class ProgressPrinter(ProgressWatcher):
	"""
		Continuously writes out the number of live and rejected points. 
	"""
	def run(self):
		import time
		while self.running:
			time.sleep(self.interval_ms / 1000.)
			if not self.running:
				break
			try:
				print(('rejected points: ', len(open(self.rejected, 'r').readlines())))
				print(('alive points: ', len(open(self.live, 'r').readlines())))
			except Exception as e:
				print(e)

class ProgressPlotter(ProgressWatcher):
	"""
		Continuously creates plots (pdfs) of the live points. 
	"""
	
	def run(self):
		import time
		while self.running:
			time.sleep(self.interval_ms / 1000.)
			if not self.running:
				break
			try:
				self._plot_live()
			except Exception as e:
				print(e)
	
	def _plot_live(self):
		import matplotlib.pyplot as plt
		import numpy
		import shutil, os
		x = numpy.loadtxt(self.live, ndmin=2)
		for i in range(self.n_params):
			plt.subplot(self.n_params, 1, 1+i)
			plt.plot(x[:,i], x[:,self.n_params], '.')
		f = "%s.pdf" % self.live
		# using a temporary file because the writing 
		# takes some time, during which the user should still be able
		# to look at the previous version
		ftmp = "%s.next.pdf" % self.live
		plt.savefig(ftmp)
		shutil.copyfile(ftmp, f)
		os.remove(ftmp)

