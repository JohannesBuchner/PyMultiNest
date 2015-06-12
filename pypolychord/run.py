from __future__ import absolute_import, unicode_literals, print_function
from ctypes import cdll
import numpy
import os

try:
	lib = cdll.LoadLibrary('libchord.so')
except OSError as e:
	if e.message == 'libchord.so: cannot open shared object file: No such file or directory':
		print()
		print('ERROR:   Could not load PolyChord library "libchord.so"')
		print('ERROR:   You have to build it first,,')
		print('ERROR:   and point the LD_LIBRARY_PATH environment variable to it!')
		print()
	if e.message.endswith('cannot open shared object file: No such file or directory'):
		print()
		print('ERROR:   Could not load PolyChord library: %s' % e.message.split(':')[0])
		print('ERROR:   You have to build PolyChord,')
		print('ERROR:   and point the LD_LIBRARY_PATH environment variable to it!')
		print()
	# the next if is useless because we can not catch symbol lookup errors (the executable crashes)
	# but it is still there as documentation.
	if 'symbol lookup error' in e.message and 'mpi' in e.message:
		print()
		print('ERROR:   You are trying to get MPI to run, but MPI failed to load.')
		print('ERROR:   Specifically, mpi symbols are missing in the executable.')
		print('ERROR:   Let me know if this is a problem of running python or a compilation problem.')
		print()
	# what if built with MPI, but don't have MPI
	print('problem:', e)
	import sys
	sys.exit(1)

from ctypes import *
import ctypes

def run(LogLikelihood, Prior, n_dims, n_live=500, n_chords=1, output_basename="chains/1-"):
	"""
	Runs PolyChord

        @param LogLikelihood:   Log-Likelihood function definition
        @param n_dims:          Number of dimensions
        @param n_live:          Number of live points
        @param n_chords:        Number of chords
        @param output_basename: File root, prefix for all output files
	
	"""

        loglike_type = CFUNCTYPE(c_double, c_int, POINTER(c_double), c_int, \
                POINTER(c_double))

        def loglike(n_dims_, theta, nderived_, phi):
            #print('py:loglike called...', n_dims_, nderived_)
            theta_pointer = cast(theta, POINTER(c_double))
            phi_pointer = cast(phi, POINTER(c_double))
            
            cube = numpy.asarray([theta_pointer[i] for i in range(n_dims)])
            params = Prior(cube)
            for i in range(n_dims):
                theta_pointer[i] = params[i]
            
            #print('py:calling likelihood...', cube, params)
            r = LogLikelihood(params)
            #print('py:returned:', r)
            return r

	"""
        loglike_type = CFUNCTYPE(c_double, POINTER(c_double), \
                POINTER(c_double), c_void_p)

        def loglike(theta, phi, context):
            theta_pointer = cast(theta, POINTER(c_double * n_dims))
            phi_pointer = cast(phi, POINTER(c_double * n_dims))
            params = np.frombuffer(theta_pointer.contents, count=n_dims)
            phi_arr = np.frombuffer(phi_pointer.contents, count=n_dims)
            
            args = [n_dims, params, phi_arr]
            
            print('calling likelihood...', args)
            r = LogLikelihood(*args)
            print('returned:', r)
            return r
        """
        
        c_double_p = ctypes.POINTER(ctypes.c_double)
        
        inifilename = output_basename + 'params.ini'
        with open(inifilename, 'w') as f:
            f.write("""
[ algorithm settings ] 
nlive = %(n_live)d
num_repeats = %(num_repeats)d
do_clustering = T
grade_frac= 1

[ posterior settings ]
weighted_posteriors = F
equally_weighted_posteriors = T
posterior_clustering = T
update_posterior = -1
boost_posterior = 5.0
base_directory = %(base)s
rootname = %(root)s

write_resume = F
resume = T
write_live = T
update_resume = %(update)d
feedback = 1

[ prior settings ]
# : name | latex name  |speed| prior type  |prior block| prior params
#--------------------------------------------------------------------
""" % (dict(n_live=n_live, 
	num_repeats=n_chords,
	base=os.path.dirname(output_basename),
	root=os.path.basename(output_basename),
	update=n_live,
	)))
	    for i in range(n_dims):
	       f.write("P : param%d | \theta_{%d} | 1 | uniform | 1 | 0 1\n" % (i, i))
        
        
        
        
        inifile = inifilename + '\0' * (100 - len(inifilename))
        Froot = output_basename + '\0' * (100 - len(output_basename))
        lib.__initsampler_MOD_dosamplingfromc(loglike_type(loglike),
                byref(c_int(n_dims)),
                create_string_buffer(inifile.encode(),100),
                None)

