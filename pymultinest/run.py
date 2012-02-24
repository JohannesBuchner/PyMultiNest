
from ctypes import cdll
lib = cdll.LoadLibrary('libcnest.so')

from ctypes import *

def run(LogLikelihood,
	Prior,
	n_dims, 
	n_params = None, 
	n_clustering_params = None, wrapped_params = None, 
	multimodal = True, const_efficiency_mode = False, n_live_points = 1000,
	evidence_tolerance = 0.5, sampling_efficiency = 0.8, 
	n_iter_before_update = 100, null_log_evidence = -1e90,
	max_modes = 100,
	outputfiles_basename = "chains/1-", seed = -1, verbose = False,
	resume = True, context = 0, write_output = True, log_zero = -1e100, 
	init_MPI = True, dump_callback = None):
	"""
	Runs MultiNest
	
	The most important parameters are the two log-probability functions Prior 
	and LogLikelihood. They are called by MultiNest.
	
	Prior should transform the unit cube into the parameter cube. Here
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
	
	Some of the parameters are explained below. Otherwise consult the 
	MultiNest documentation.
	
	@param n_params: 
		Total no. of parameters, should be equal to ndims in most cases 
		but if you need to store some additional
		parameters with the actual parameters then you need to pass 
		them through the likelihood routine.

	@param sampling_efficiency:
		defines the sampling efficiency. 0.8 and 0.3 are recommended 
		for parameter estimation & evidence evalutation
		respectively.
		use 'parameter' or 'model' to select the respective default 
		values


	@param evidence_tolerance:
		A value of 0.5 should give good enough accuracy.

	@param n_clustering_params:
		If mmodal is T, MultiNest will attempt to separate out the 
		modes. Mode separation is done through a clustering
		algorithm. Mode separation can be done on all the parameters 
		(in which case nCdims should be set to ndims) & it
		can also be done on a subset of parameters (in which case 
		nCdims < ndims) which might be advantageous as
		clustering is less accurate as the dimensionality increases. 
		If nCdims < ndims then mode separation is done on
		the first nCdims parameters.

	@param null_log_evidence:
		If mmodal is T, MultiNest can find multiple modes & also specify 
		which samples belong to which mode. It might be
		desirable to have separate samples & mode statistics for modes 
		with local log-evidence value greater than a
		particular value in which case nullZ should be set to that 
		value. If there isn't any particulrly interesting
		nullZ value, then nullZ should be set to a very large negative 
		number (e.g. -1.d90).
		
	@param init_MPI:
		initialize MPI routines?, relevant only if compiling with MPI
	
	@param log_zero: 
		points with loglike < logZero will be ignored by MultiNest
	
	@param write_output:
		write output files? This is required for analysis.
		
	@param dump_callback:
		a callback function for dumping the current status
	
	"""

	if n_params == None:
		n_params = n_dims
	if n_clustering_params == None:
		n_clustering_params = n_params
	if wrapped_params == None:
		wrapped_params = [0] * n_params
	
	WrappedType = c_int * len(wrapped_params)
	wraps = WrappedType(*wrapped_params)
	
	if sampling_efficiency == 'parameter':
		sampling_efficiency = 0.8
	if sampling_efficiency == 'model':
		sampling_efficiency = 0.3
	
	
	prior_type = CFUNCTYPE(c_void_p, POINTER(c_double), c_int, c_int)
	loglike_type = CFUNCTYPE(c_double, POINTER(c_double), c_int, c_int)
	dumper_type = CFUNCTYPE(c_void_p, 
		c_int, c_int, c_int, POINTER(c_double))
	""", 
		POINTER(c_double), POINTER(c_double), POINTER(c_double), 
		POINTER(c_double), POINTER(c_double), POINTER(c_double),
		c_double, c_double, c_double)"""
	lib.reset()
	
	lib.set_function(prior_type(Prior), loglike_type(LogLikelihood))
	
	if dump_callback is not None:
		lib.set_dumper(dumper_type(dump_wrapper))

	lib.run(c_int(multimodal), c_int(const_efficiency_mode), 
		c_int(n_live_points), c_double(evidence_tolerance), 
		c_double(sampling_efficiency), c_int(n_dims), c_int(n_params),
		c_int(n_clustering_params), c_int(max_modes), 
		c_int(n_iter_before_update), c_double(evidence_tolerance), 
		outputfiles_basename, c_int(seed), wraps,
		c_int(verbose), c_int(resume), 
		c_int(write_output), c_int(init_MPI), 
		c_double(log_zero), c_int(context))


