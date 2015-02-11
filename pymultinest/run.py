from __future__ import absolute_import, unicode_literals, print_function
from ctypes import cdll
import sys

libname = 'libmultinest'
try: # detect if run through mpiexec/mpirun
	from mpi4py import MPI
	if MPI.COMM_WORLD.Get_size() > 1: # need parallel capabilities
		libname = 'libmultinest_mpi'
except ImportError:
	pass

libname += {
	'darwin' : '.dylib',
	'win32'  : '.dll',
	'cygwin' : '.dll',
}.get(sys.platform, '.so')

try:
	lib = cdll.LoadLibrary(libname)
except OSError as e:
	if e.message == '%s: cannot open shared object file: No such file or directory' % libname:
		print()
		print('ERROR:   Could not load MultiNest library "%s"' % libname)
		print('ERROR:   You have to build it first,')
		print('ERROR:   and point the LD_LIBRARY_PATH environment variable to it!')
		print('ERROR:   manual: http://johannesbuchner.github.com/PyMultiNest/install.html')
		print()
	if e.message.endswith('cannot open shared object file: No such file or directory'):
		print()
		print('ERROR:   Could not load MultiNest library: %s' % e.message.split(':')[0])
		print('ERROR:   You have to build MultiNest,')
		print('ERROR:   and point the LD_LIBRARY_PATH environment variable to it!')
		print('ERROR:   manual: http://johannesbuchner.github.com/PyMultiNest/install.html')
		print()
	if 'undefined symbol: mpi_' in e.message:
		print()
		print('ERROR:   You tried to compile MultiNest linked with MPI,')
		print('ERROR:   but now when running, MultiNest can not find the MPI linked libraries.')
		print('ERROR:   manual: http://johannesbuchner.github.com/PyMultiNest/install.html')
		print()
	# the next if is useless because we can not catch symbol lookup errors (the executable crashes)
	# but it is still there as documentation.
	if 'symbol lookup error' in e.message and 'mpi' in e.message:
		print()
		print('ERROR:   You are trying to get MPI to run, but MPI failed to load.')
		print('ERROR:   Specifically, mpi symbols are missing in the executable.')
		print('ERROR:   Let me know if this is a problem of running python or a compilation problem.')
		print('ERROR:   manual: http://johannesbuchner.github.com/PyMultiNest/install.html')
		print()
	# what if built with MPI, but don't have MPI
	print('problem:', e)
	import sys
	sys.exit(1)

from ctypes import *
from numpy.ctypeslib import as_array
import signal, sys
import inspect

def interrupt_handler(signal, frame):
	sys.stderr.write('ERROR: Interrupt received: Terminating\n')
	sys.exit(1)

def run(LogLikelihood,
	Prior,
	n_dims, 
	n_params = None, 
	n_clustering_params = None, wrapped_params = None, 
	importance_nested_sampling = True,
	multimodal = True, const_efficiency_mode = False, n_live_points = 400,
	evidence_tolerance = 0.5, sampling_efficiency = 0.8, 
	n_iter_before_update = 100, null_log_evidence = -1e90,
	max_modes = 100, mode_tolerance = -1e90,
	outputfiles_basename = "chains/1-", seed = -1, verbose = False,
	resume = True, context = 0, write_output = True, log_zero = -1e100, 
	max_iter = 0, init_MPI = False, dump_callback = None):
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
	
		def Loglike(cube, ndim, nparams, lnew):
			chi = 1.
			
			for i in range(ndim):
				chi *= math.cos(cube[i] / 2.)
			return math.pow(2. + chi, 5)
	
	Some of the parameters are explained below. Otherwise consult the 
	MultiNest documentation.
	
	@param importance_nested_sampling:
		If True, Multinest will use Importance Nested Sampling (INS). Read http://arxiv.org/abs/1306.2144
		for more details on INS. Please read the MultiNest README file before using the INS in MultiNest v3.0.
	
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
		To run pymultinest with MPI, you need mpi4py installed. Then,
		the libmultinest_mpi library is loaded when you run with mpiexec
		or similar. init_MPI should be set to False, because importing
		mpi4py initialises MPI already.
	
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

	loglike_type = CFUNCTYPE(c_double, POINTER(c_double),
		c_int, c_int, c_double, c_void_p)

	dumper_type  = CFUNCTYPE(c_void_p, c_int, c_int, c_int,
		POINTER(c_double),POINTER(c_double),POINTER(c_double),
		c_double,c_double,c_double,c_void_p)
	
	# check if lnew is supported by user function
	nargs = 3
	try:
		nargs = len(inspect.getargspec(LogLikelihood).args) - inspect.ismethod(LogLikelihood)
	except:
		pass
	
	if nargs == 4:
		def loglike(cube, ndim, nparams, lnew, nullcontext):
			if Prior:
				Prior(cube, ndim, nparams)
			return LogLikelihood(cube, ndim, nparams, lnew)
	else:
		def loglike(cube, ndim, nparams, lnew, nullcontext):
			if Prior:
				Prior(cube, ndim, nparams)
			return LogLikelihood(cube, ndim, nparams)
	
	def dumper(nSamples,nlive,nPar,
			   physLive,posterior,paramConstr,
			   maxLogLike,logZ,logZerr,nullcontext):
		if dump_callback:
			# It's not clear to me what the desired PyMultiNest dumper callback
			# syntax is... but this should pass back the right numpy arrays,
			# without copies. Untested!
			pc =  as_array(paramConstr,shape=(nPar,4))
			
			dump_callback(nSamples,nlive,nPar,
				as_array(physLive,shape=(nPar+1,nlive)).T,
				as_array(posterior,shape=(nPar+2,nSamples)).T,
				(pc[0,:],pc[1,:],pc[2,:],pc[3,:]), # (mean,std,bestfit,map)
				maxLogLike,logZ,logZerr)
	prev_handler = signal.signal(signal.SIGINT, interrupt_handler)
	
	lib.run(c_bool(importance_nested_sampling),
		c_bool(multimodal), c_bool(const_efficiency_mode),
		c_int(n_live_points), c_double(evidence_tolerance), 
		c_double(sampling_efficiency), c_int(n_dims), c_int(n_params),
		c_int(n_clustering_params), c_int(max_modes), 
		c_int(n_iter_before_update), c_double(mode_tolerance), 
		create_string_buffer(outputfiles_basename.encode(),100),
		c_int(seed), wraps,
		c_bool(verbose), c_bool(resume),
		c_bool(write_output), c_bool(init_MPI),
		c_double(log_zero), c_int(max_iter),
		loglike_type(loglike),dumper_type(dumper),
		c_int(context))
	signal.signal(signal.SIGINT, prev_handler)


