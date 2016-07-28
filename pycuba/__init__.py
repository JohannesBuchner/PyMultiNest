from __future__ import absolute_import, unicode_literals, print_function
import os
import ctypes
from ctypes import POINTER, c_int, c_double, c_void_p, byref

"""
Parallelisation within cuba is not supported, because python does not know
that the call is in parallel and writes to the same memory location, causing
overrides. This could be overcome by using locks.
"""
os.environ['CUBACORES'] = '0'
lib = ctypes.cdll.LoadLibrary('libcuba.so')

NULL = ctypes.POINTER(c_int)()

# defaults
EPSREL = 1e-3
EPSABS = 1e-12
LAST = 4
MINEVAL = 0
MAXEVAL = 50000

NSTART = 1000
NINCREASE = 500
NBATCH = 1000
GRIDNO = 0
STATEFILE = NULL
spin = NULL

class BOUNDS(ctypes.Structure):
  _fields_ = ("lower", c_double), ("upper", c_double)

integrand_type = ctypes.CFUNCTYPE(c_int, POINTER(c_int), 
  POINTER(c_double), POINTER(c_int), POINTER(c_double), c_void_p)
peakfinder_type = ctypes.CFUNCTYPE(c_void_p, POINTER(c_int), POINTER(BOUNDS), 
  POINTER(c_int), POINTER(c_double))

def wrap_integrand(integrand):
  use_raw_callback = isinstance(integrand, ctypes._CFuncPtr)
  if use_raw_callback:
    return integrand
  else:
    return integrand_type(integrand)
  
def Vegas(integrand, ndim, userdata=NULL, 
    epsrel=EPSREL, epsabs=EPSABS, verbose=0, ncomp=1, seed=None,
    mineval=MINEVAL, maxeval=MAXEVAL, nstart=NSTART, 
    nincrease=NINCREASE, nbatch=NBATCH,
    gridno=GRIDNO, statefile=NULL, nvec=1):
  """
  *nstart*: the number of integrand evaluations per iteration to start
with.
  *nincrease*: the increase in the number of integrand evaluations per
iteration.
  *nbatch*: the batch size for sampling.
    Vegas samples points not all at once, but in batches of size nbatch, to avoid exces-
    sive memory consumption. 1000 is a reasonable value, though it should not affect
    performance too much.
  *gridno*: the slot in the internal grid table.
    It may accelerate convergence to keep the grid accumulated during one integration for
    the next one, if the integrands are reasonably similar to each other. Vegas maintains
    an internal table with space for ten grids for this purpose. The slot in this grid is
    specified by gridno.
    If a grid number between 1 and 10 is selected, the grid is not discarded at the end of
    the integration, but stored in the respective slot of the table for a future invocation.
    The grid is only re-used if the dimension of the subsequent integration is the same
    as the one it originates from.
    In repeated invocations it may become necessary to flush a slot in memory, in which
    case the negative of the grid number should be set.
  *statefile*: a filename for storing the internal state.
    Vegas can store its entire internal state (i.e. all the information to resume an inter-
    rupted integration) in an external file. The state file is updated after every iteration.
    If, on a subsequent invocation, Vegas finds a file of the specified name, it loads the
    internal state and continues from the point it left off. Needless to say, using an ex-
    isting state file with a different integrand generally leads to wrong results. Once the
    integration finishes successfully, i.e. the prescribed accuracy is attained, the state file
    is removed.
    This feature is useful mainly to define 'check-points' in long-running integrations
    from which the calculation can be restarted.
  
  Vegas actually passes the integrand two more arguments, i.e. the integrand subroutine
  is really declared as
  
  def integrand(ndim, x, ncomp, f, userdata, weight, iter):
  
  where weight contains the weight of the point being sampled and iter the current iteration
  number. These extra arguments can safely be ignored, i.e. the integrand may be (and
  usually is) defined with just five (or even four, if userdata is also ignored) arguments.
  """
  
  neval = c_int()
  fail = c_int()
  comp = c_int()
  ARR = c_double * ncomp
  integral = ARR()
  error = ARR()
  prob  = ARR()
  
  if seed is None:
    seed = 0
  
  lib.Vegas(ndim, ncomp, wrap_integrand(integrand), userdata,
    c_int(nvec), c_double(epsrel), c_double(epsabs), verbose, seed,
    mineval, maxeval, nstart, nincrease, nbatch,
    gridno, statefile, spin,
    byref(neval), byref(fail), integral, error, prob)
  
  return dict(neval=neval.value, fail=fail.value, comp=comp.value,
    results=[{
      'integral':integral[comp], 
      'error':error[comp], 
      'prob':prob[comp]
      } for comp in range(ncomp)])

def Suave(integrand, ndim, nnew=1000, nmin=2, flatness=50., userdata=NULL, 
    epsrel=EPSREL, epsabs=EPSABS, verbose=0, ncomp=1, seed=None,
    mineval=MINEVAL, maxeval=MAXEVAL, statefile=NULL, nvec=1):
  """
  *nnew*: the number of new integrand evaluations in each subdivision.
  
  *nmin*: the minimum number of samples a former pass must contribute
  to a subregion to be considered in that region's compound integral value. Increasing
  nmin may reduce jumps in the chi^2 value.
  
  *flatness*: the parameter p in Eq. (1), i.e. the type of norm
  used to compute the fluctuation of a sample. This determines how prominently 'out-
  liers,' i.e. individual samples with a large fluctuation, figure in the total fluctuation,
  which in turn determines how a region is split up. As suggested by its name, flatness
  should be chosen large for 'flat' integrands and small for 'volatile' integrands with
  high peaks. Note that since flatness appears in the exponent, one should not use
  too large values (say, no more than a few hundred) lest terms be truncated internally
  to prevent overflow.
  """
  
  neval = c_int()
  fail = c_int()
  comp = c_int()
  nregions = c_int()
  ARR = c_double * ncomp
  integral = ARR()
  error = ARR()
  prob  = ARR()
    
  if seed is None:
    seed = 0
  
  lib.Suave(ndim, ncomp, wrap_integrand(integrand), userdata,
    c_int(nvec), c_double(epsrel), c_double(epsabs), verbose, seed,
    mineval, maxeval, nnew, nmin, c_double(flatness), statefile, spin,
    byref(nregions), byref(neval), byref(fail), integral, error, prob)
  
  return dict(neval=neval.value, fail=fail.value, comp=comp.value, nregions=nregions.value,
    results=[{
      'integral':integral[comp], 
      'error':error[comp], 
      'prob':prob[comp]
      } for comp in range(ncomp)])
  
def Divonne(integrand, ndim, 
    key1, key2, key3, maxpass, border,
    maxchisq, mindeviation,
    mineval=MINEVAL, maxeval=MAXEVAL, ncomp=1,
    ldxgiven=None, xgiven=None, nextra=0, peakfinder=None,
    userdata=NULL, seed=None,
    epsrel=EPSREL, epsabs=EPSABS, verbose=0, statefile=NULL, nvec=1):
  """
  *key1*: determines sampling in the partitioning phase:
    key1 = 7, 9, 11, 13 selects the cubature rule of degree key1. Note that the degree-11
    rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
    For other values of key1, a quasi-random sample of n1 = |key1| points is used, where
    the sign of key1 determines the type of sample,
    - key1 > 0, use a Korobov quasi-random sample,
    - key1 < 0, use a "standard" sample (a Sobol quasi-random sample if seed = 0,
    otherwise a pseudo-random sample).
  
  *key2*: determines sampling in the final integration phase:
    key2 = 7, 9, 11, 13 selects the cubature rule of degree key2. Note that the degree-11
    rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
    For other values of key2, a quasi-random sample is used, where the sign of key2
    determines the type of sample,
    - key2 > 0, use a Korobov quasi-random sample,
    - key2 < 0, use a 'standard' sample (see description of key1 above),
    and n2 = |key2| determines the number of points,
    - n2 >= 40, sample n2 points,
    - n2 < 40, sample n2 need points, where nneed is the number of points needed to
    reach the prescribed accuracy, as estimated by Divonne from the results of the
    partitioning phase.
    
  *key3* sets the strategy for the refinement phase:
    key3 = 0, do not treat the subregion any further.
    key3 = 1, split the subregion up once more.
    Otherwise, the subregion is sampled a third time with key3 specifying the sampling
    parameters exactly as key2 above.
  
  *maxpass*: controls the thoroughness of the partitioning phase: The
    partitioning phase terminates when the estimated total number of integrand evalu-
    ations (partitioning plus final integration) does not decrease for maxpass successive
    iterations.
    A decrease in points generally indicates that Divonne discovered new structures of
    the integrand and was able to find a more effective partitioning. maxpass can be
    understood as the number of 'safety' iterations that are performed before the par-
    tition is accepted as final and counting consequently restarts at zero whenever new
    structures are found.
    
  *border*: the width of the border of the integration region.
    Points falling into this border region will not be sampled directly, but will be extrap-
    olated from two samples from the interior. Use a non-zero border if the integrand
    subroutine cannot produce values directly on the integration boundary.
  
  *maxchisq*: the maximum chisquare value a single subregion is al-
    lowed to have in the final integration phase. Regions which fail this chisquare test and whose
    sample averages differ by more than mindeviation move on to the refinement phase.
 
  *mindeviation*: a bound, given as the fraction of the re-
    quested error of the entire integral, which determines whether it is worthwhile fur-
    ther examining a region that failed the chisquare test. Only if the two sampling averages
    obtained for the region differ by more than this bound is the region further treated.
 
 *ldxgiven*: the leading dimension of xgiven, i.e. the offset between one
    point and the next in memory.
 
 *xgiven(ldxgiven,ngiven)*: a list of points where the inte-
    grand might have peaks. Divonne will consider these points when partitioning the
    integration region. The idea here is to help the integrator find the extrema of the in-
    tegrand in the presence of very narrow peaks. Even if only the approximate location
    of such peaks is known, this can considerably speed up convergence.
 
 *nextra*: the maximum number of extra points the peak-finder subrou-
    tine will return. If nextra is zero, peakfinder is not called and an arbitrary object
    may be passed in its place, e.g. just 0.
  
 *peakfinder*: the peak-finder subroutine. This subroutine is called
    whenever a region is up for subdivision and is supposed to point out possible peaks
    lying in the region, thus acting as the dynamic counterpart of the static list of points
    supplied in xgiven. It is expected to be declared as
    def peakfinder(ndim, b, n, x):
      The bounds of the subregion are passed in the array b, where b(1,d ) is the lower and
      b(2,d ) the upper bound in dimension d . On entry, n specifies the maximum number
      of points that may be written to x. On exit, n must contain the actual number of
      points in x.

  Divonne actually passes the integrand one more argument, i.e. the integrand subroutine is
  really declared as
  
  def integrand(ndim, x, ncomp, f, phase):
  
  The fifth argument, phase, indicates the integration phase:
   0, sampling of the points in xgiven,
   1, partitioning phase,
   2, final integration phase,
   3, refinement phase.
   
  This information might be useful if the integrand takes long to compute and a sufficiently
  accurate approximation of the integrand is available. The actual value of the integral is only
  of minor importance in the partitioning phase, which is instead much more dependent on
  the peak structure of the integrand to find an appropriate tessellation. An approximation
  which reproduces the peak structure while leaving out the fine details might hence be a
  perfectly viable and much faster substitute when phase < 2.
  
  In all other instances, phase can be ignored and it is entirely admissible to define the
  integrand with only five arguments.
  """
  
  neval = c_int()
  fail = c_int()
  comp = c_int()
  nregions = c_int()
  ARR = c_double * ncomp
  integral = ARR()
  error = ARR()
  prob  = ARR()
  if ldxgiven is None:
    ldxgiven = ndim
  
  if xgiven is not None:
    ngiven = len(xgiven)
    xgiven = ARR(xgiven)
  else:
    ngiven = 0
    xgiven = NULL
  
  if peakfinder is None:
    peakfinder = NULL
  else:
    peakfinder = peakfinder_type(peakfinder)
  
  if seed is None:
    seed = 0

  lib.Divonne(ndim, ncomp, wrap_integrand(integrand), userdata,
    c_int(nvec), c_double(epsrel), c_double(epsabs), verbose, seed,
    mineval, maxeval, key1, key2, key3, maxpass, 
    c_double(border), c_double(maxchisq), c_double(mindeviation), 
    ngiven, ldxgiven, xgiven, nextra, peakfinder, statefile, spin,
    byref(nregions), byref(neval), byref(fail), integral, error, prob)
  
  return dict(neval=neval.value, fail=fail.value, comp=comp.value, nregions=nregions.value,
    results=[{
      'integral':integral[comp], 
      'error':error[comp], 
      'prob':prob[comp]
      } for comp in range(ncomp)])

def Cuhre(integrand, ndim, 
    key=0, mineval=MINEVAL, maxeval=MAXEVAL, ncomp=1,
    userdata=NULL, seed=None,
    epsrel=EPSREL, epsabs=EPSABS, verbose=0, statefile=NULL, nvec=1):
  """
  *key* chooses the basic integration rule:
    key = 7, 9, 11, 13 selects the cubature rule of degree key. Note that the degree-11
    rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
    For other values, the default rule is taken, which is the degree-13 rule in 2 dimensions,
    the degree-11 rule in 3 dimensions, and the degree-9 rule otherwise.
  """
  
  neval = c_int()
  fail = c_int()
  comp = c_int()
  nregions = c_int()
  ARR = c_double * ncomp
  integral = ARR()
  error = ARR()
  prob  = ARR()
  
  if seed is None:
    seed = 0

  lib.Cuhre(ndim, ncomp, wrap_integrand(integrand), userdata,
    c_int(nvec), c_double(epsrel), c_double(epsabs), verbose,
    mineval, maxeval, key, statefile, spin,
    byref(nregions), byref(neval), byref(fail), integral, error, prob)
  
  return dict(neval=neval.value, fail=fail.value, comp=comp.value, nregions=nregions.value,
    results=[{
      'integral':integral[comp], 
      'error':error[comp], 
      'prob':prob[comp]
      } for comp in range(ncomp)])

def demo():
  import math

  def Integrand(ndim, xx, ncomp, ff, userdata):
    x,y,z = [xx[i] for i in range(ndim.contents.value)]
    result = math.sin(x)*math.cos(y)*math.exp(z)
    ff[0] = result
    return 0
    

  NDIM = 3
  NCOMP = 1

  NNEW = 1000
  NMIN = 2
  FLATNESS = 50.

  KEY1 = 47
  KEY2 = 1
  KEY3 = 1
  MAXPASS = 5
  BORDER = 0.
  MAXCHISQ = 10.
  MINDEVIATION = .25
  NGIVEN = 0
  LDXGIVEN = NDIM
  NEXTRA = 0
  MINEVAL = 0
  MAXEVAL = 50000


  KEY = 0
  
  def print_header(name):
    print('-------------------- %s test -------------------' % name)
  def print_results(name, results):
    keys = ['nregions', 'neval', 'fail']
    keys = list(filter(results.has_key, keys))
    text = ["%s %d" % (k, results[k]) for k in keys]
    print("%s RESULT:\t" % name.upper() + "\t".join(text))
    for comp in results['results']:
      print("%s RESULT:\t" % name.upper() + \
	"%(integral).8f +- %(error).8f\tp = %(prob).3f\n" % comp)
  
  from os import environ as env
  verbose = 2
  if 'CUBAVERBOSE' in env:
    verbose = int(env['CUBAVERBOSE'])
  
  print_header('Vegas')
  print_results('Vegas', Vegas(Integrand, NDIM, verbose=2))
  
  print_header('Suave')
  print_results('Suave', Suave(Integrand, NDIM, NNEW, FLATNESS, verbose=2 | 4))
  
  print_header('Divonne')
  print_results('Divonne', Divonne(Integrand, NDIM, 
    mineval=MINEVAL, maxeval=MAXEVAL,
    key1=KEY1, key2=KEY2, key3=KEY3, maxpass=MAXPASS,
    border=BORDER, maxchisq=MAXCHISQ, mindeviation=MINDEVIATION,
    ldxgiven=LDXGIVEN, verbose=2))
  
  print_header('Cuhre')
  print_results('Cuhre', Cuhre(Integrand, NDIM, key=KEY, verbose=2 | 4))
  
if __name__ == '__main__':
  demo()



