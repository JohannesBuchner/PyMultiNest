"""
Module for watching the progress of APEMoST
"""

import threading
import analyse

class ProgressWatcher(threading.Thread):
	"""
		Watches the progress of APEMoST.
		
		Not implemented yet.
	"""
	def __init__(self, interval_ms = 200, outputfiles_basename = "."):
		threading.Thread.__init__(self)
		self.params = analyse.load_params()
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
				import traceback
				traceback.print_exc()
	
	def stop(self):
		self.running = False
	
	def _plot_live(self):
		# there are some issues because simultaneous read access to the output files 
		# seems to cause segmentation faults.
		return
	
