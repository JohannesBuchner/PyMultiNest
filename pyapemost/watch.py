"""
Module for watching the progress of APEMoST
"""
from __future__ import absolute_import, unicode_literals, print_function
import threading
from . import analyse
import matplotlib.pyplot as plt
import shutil, os

class ProgressWatcher(threading.Thread):
	"""
		Watches the progress of APEMoST.
		
		Not implemented yet.
	"""
	def __init__(self, interval_ms = 200, outputfiles_basename = "", window = False):
		threading.Thread.__init__(self)
		self.params = analyse.load_params()
		self.outputfiles_basename = outputfiles_basename
		self.interval_ms = interval_ms
		self.running = True
		
		if window:
			self.plotter = analyse.VisitedWindow()
			self.tmpfile = False
		else:
			self.plotter = analyse.VisitedAllPlotter(self.outputfiles_basename + "live-next-", nlast = 1000)
			self.tmpfile = self.plotter.outputfiles_basename + "chain0.pdf"
			self.outputfile = self.outputfiles_basename + "live-chain0.pdf"

	def run(self):
		import time
		while self.running:
			time.sleep(self.interval_ms / 1000.)
			if not self.running:
				break
			try:
				self._plot_live()
			except Exception as e:
				print('Plotting failed', e)
				#import traceback
				#traceback.print_exc()
	
	def stop(self):
		self.running = False
	
	def _plot_live(self):
		self.plotter.plot()
		if self.outputfile:
			# using a temporary file because the writing 
			# takes some time, during which the user should still be able
			# to look at the previous version
			shutil.copyfile(self.tmpfile, self.outputfile)
			os.remove(self.tmpfile)
		
	
