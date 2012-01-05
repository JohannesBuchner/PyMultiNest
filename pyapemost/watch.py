"""
Module for watching the progress of APEMoST
"""

import threading

class ProgressWatcher(threading.Thread):
	"""
		Watches the progress of APEMost.
		
		Not implemented yet.
	"""
	def __init__(self, n_params, interval_ms = 200, outputfiles_basename = "."):
		threading.Thread.__init__(self)
		self.n_params = n_params
		self.outputfiles_basename = outputfiles_basename
		self.interval_ms = interval_ms
		self.running = True

	def run(self):
		import time
		while self.running:
			time.sleep(self.interval_ms / 1000.)
			if not self.running:
				break
			try:
				self._plot_live()
			except Exception as e:
				print e
	
	def _plot_live(self):
		pass
	

