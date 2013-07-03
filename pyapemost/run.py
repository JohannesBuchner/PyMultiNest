from __future__ import absolute_import, unicode_literals, print_function
from ctypes import *
try:
	lib = cdll.LoadLibrary('libapemost.so')
except OSError as e:
	if e.message == 'libapemost.so: cannot open shared object file: No such file or directory':
		print()
		print('ERROR:   Could not load apemost library "libapemost.so"')
		print('ERROR:   You have to build it (download APEMoST), and point ')
		print('ERROR:   the LD_LIBRARY_PATH environment variable to it!')
		print()
	print(e)
	import sys
	sys.exit(1)

import ctypesGsl as gsl

"""
	The current state of the Markov Chain. Foremost, you need to
	access params using ctypesGsl.vector functions.
	Secondly, data contains the file "data" read in as a matrix.
"""
class MCMC(Structure):
	_fields_ = [
		('n_params', c_uint),
		('accept', c_ulong),
		('reject', c_ulong),
		('probability', c_double),
		('prior', c_double),
		('probability_best', c_double),
		('rng', POINTER(gsl.gsl_rng)),
		('params', POINTER(gsl.gsl_vector)),
		('params_best', POINTER(gsl.gsl_vector)),
		('dice', c_double),
		('files', POINTER(c_void_p)),
		('parameter_names', POINTER(c_char_p)),
		('accept', POINTER(c_ulong)),
		('reject', POINTER(c_ulong)),
		('params_stepwidth', POINTER(gsl.gsl_vector)),
		('params_min', POINTER(gsl.gsl_vector)),
		('params_max', POINTER(gsl.gsl_vector)),
		('data', POINTER(gsl.gsl_matrix)),
		('n_iterations', c_ulong),
		('extra_data', c_void_p),
	]
	
	def get_params():
		return (self.params.contents.data[i] for i in range(m.n_params))


callback_type = CFUNCTYPE(c_double, POINTER(MCMC), POINTER(gsl.gsl_vector))
_last_functions = None
def set_function(LogLikelihood, Prior):
	"""
	Initializes APEMoST with the two log-probability functions 
	Prior and LogLikelihood. They are called by APEMoST.
	
	Prior is given the current parameters and should return a
	log
	should transform the unit cube into the parameter cube. Here
	is an example for a uniform prior::
	
		def Prior(cube, ndim, nparams):
			for i in range(ndim):
				cube[i] = cube[i] * 10 * math.pi
	
	The LogLikelihood function gets this parameter cube and should
	return the logarithm of the likelihood.
	Here is the example for the eggbox problem::
	
		def Loglike(cube, ndim, nparams):
			chi = 1.
			
			for i in range(ndim):
				chi *= math.cos(cube[i] / 2.)
			return math.pow(2. + chi, 5)
	@param LogLikelihood: Log of the Likelihood function
		(without Prior)

	@param Prior: Log of the Prior function
	"""
	global _last_functions
	# avoid garbage collection of functions, which leads to segfaults
	_last_functions = callback_type(LogLikelihood), callback_type(Prior)
	lib.set_function(*_last_functions)
	

calibrate_first_chain = lib.calibrate_first
"""
Calibrate the first chain
"""
calibrate_other_chains = lib.calibrate_rest
"""
Calibrate the other chains
"""

def calibrate():
	"""
	Calibrate all chains
	"""
	calibrate_first_chain()
	calibrate_other_chains()

	
def run(max_iterations = 0, append = True):
	"""
	Run the MCMC sampler.
	
	@param max_iterations: stop after this many iterations.
	if 0, run indefinitely (until SIGTERM received).
	
	@param append: Append to existing files if found?
	"""
	lib.prepare_and_run_sampler(max_iterations, append)



