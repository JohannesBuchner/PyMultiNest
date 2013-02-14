
from ctypes import cdll
try:
	lib = cdll.LoadLibrary('libcnest.so')
except OSError as e:
	if e.message == 'libcnest.so: cannot open shared object file: No such file or directory':
		print
		print 'ERROR:   Could not load MultiNest Bridge library "libcnest.so"'
		print 'ERROR:   You have to build it (from the multinest_bridge folder in pymultinest),'
		print 'ERROR:   and point the LD_LIBRARY_PATH environment variable to it!'
		print 'ERROR:   manual: http://johannesbuchner.github.com/PyMultiNest/install.html'
		print
	if e.message == 'libnest3.so: cannot open shared object file: No such file or directory':
		print
		print 'ERROR:   Could not load MultiNest library "libnest3.so"'
		print 'ERROR:   You have to build it (in MultiNest, run make libnest3.so WITHOUT_MPI=1),'
		print 'ERROR:   and point the LD_LIBRARY_PATH environment variable to it!'
		print 'ERROR:   manual: http://johannesbuchner.github.com/PyMultiNest/install.html'
		print
	if 'undefined symbol: mpi_' in e.message:
		print
		print 'ERROR:   You did something stupid. You tried to compile MultiNest with MPI,'
		print 'ERROR:   but the MultiNest bridge without MPI, or the other way around.'
		print 'ERROR:   Decide on one consistent way (no MPI is safe). Either way, you currently'
		print 'ERROR:   do not have MPI compiled into the library, and it fails to load!'
		print 'ERROR:   manual: http://johannesbuchner.github.com/PyMultiNest/install.html'
		print
	if 'libcnest.so: undefined symbol:' in e.message and 'nestrun' in e.message:
		print
		print 'ERROR:   Sorry.'
		print 'ERROR:   The compiler decided to call the MultiNest routine differently than we expect.'
		print 'ERROR:   You have to run $ readelf -s libnest3.so | grep nestrun'
		print 'ERROR:   to find out how the nestrun routine is called (typically __nested_MOD_nestrun).'
		print 'ERROR:   Then you have to compile the MultiNest bridge using'
		print 'ERROR:   $ make clean libcnest.so MULTINEST_CALL=__nested_MOD_nestrun '
		print 'ERROR:   manual: http://johannesbuchner.github.com/PyMultiNest/install.html'
		print
	# the next if is useless because we can not catch symbol lookup errors (the executable crashes)
	# but it is still there as documentation.
	if 'symbol lookup error' in e.message and 'mpi' in e.message:
		print
		print 'ERROR:   You are trying to get MPI to run, but MPI failed to load.'
		print 'ERROR:   Specifically, mpi symbols are missing in the executable.'
		print 'ERROR:   Let me know if this is a problem of running python or a compilation problem.'
		print 'ERROR:   manual: http://johannesbuchner.github.com/PyMultiNest/install.html'
		print
	# what if built with MPI, but don't have MPI
	print e
	import sys
	sys.exit(1)

from ctypes import *

def run(LogLikelihood,
	Prior,
	n_dims, 
	n_params = None, 
	n_clustering_params = None, wrapped_params = None, 
	multimodal = True, const_efficiency_mode = False, n_live_points = 1000,
	evidence_tolerance = 0.5, sampling_efficiency = 0.8, 
	n_iter_before_update = 100, null_log_evidence = -1e90,
	max_modes = 100, mode_tolerance = 1e-90,
	outputfiles_basename = "chains/1-", seed = -1, verbose = False,
	resume = True, context = 0, write_output = True, log_zero = -1e100, 
	max_iter = 0, init_MPI = True, dump_callback = None):
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

	@param mode_tolerance:
		MultiNest can find multiple modes & also specify which samples belong to which mode. It might be
		desirable to have separate samples & mode statistics for modes with local log-evidence value greater than a
		particular value in which case Ztol should be set to that value. If there isn't any particularly interesting
		Ztol value, then Ztol should be set to a very large negative number (e.g. -1e90).

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
	
	@param max_iter: 
		maximum number of iterations. 0 is unlimited.
	
	@param write_output:
		write output files? This is required for analysis.
		
	@param dump_callback:
		a callback function for dumping the current status
	
	"""

	if n_params == None:
		n_params = n_dims
	if n_clustering_params == None:
		n_clustering_params = n_dims
	if wrapped_params == None:
		wrapped_params = [0] * n_dims
	
	WrappedType = c_int * len(wrapped_params)
	wraps = WrappedType(*wrapped_params)
	
	if sampling_efficiency == 'parameter':
		sampling_efficiency = 0.8
	if sampling_efficiency == 'model':
		sampling_efficiency = 0.3
	
	lib.reset()
	
	prior_type = CFUNCTYPE(c_void_p, POINTER(c_double), c_int, c_int)
	if Prior is not None:
		c_Prior = prior_type(Prior)
		lib.set_prior(c_Prior)
	
	loglike_type = CFUNCTYPE(c_double, POINTER(c_double), c_int, c_int)
	c_LogLikelihood = loglike_type(LogLikelihood)
	lib.set_function(c_LogLikelihood)
	
	dumper_type = CFUNCTYPE(c_void_p, 
		c_int, c_int, c_int, POINTER(c_double))
	""", 
		POINTER(c_double), POINTER(c_double), POINTER(c_double), 
		POINTER(c_double), POINTER(c_double), POINTER(c_double),
		c_double, c_double, c_double)"""
	
	if dump_callback is not None:
		lib.set_dumper(dumper_type(dump_wrapper))

	lib.run(c_int(multimodal), c_int(const_efficiency_mode), 
		c_int(n_live_points), c_double(evidence_tolerance), 
		c_double(sampling_efficiency), c_int(n_dims), c_int(n_params),
		c_int(n_clustering_params), c_int(max_modes), 
		c_int(n_iter_before_update), c_double(mode_tolerance), 
		outputfiles_basename, c_int(seed), wraps,
		c_int(verbose), c_int(resume), 
		c_int(write_output), c_int(init_MPI), 
		c_double(log_zero), c_int(max_iter),
		c_int(context))


