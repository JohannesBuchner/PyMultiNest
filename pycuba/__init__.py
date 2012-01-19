import ctypes
from ctypes import POINTER, c_int, c_double, c_void_p, byref

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

class BOUNDS(ctypes.Structure):
  _fields_ = ("lower", c_double), ("upper", c_double)

integrand_type = ctypes.CFUNCTYPE(c_int, POINTER(c_int), 
  POINTER(c_double), POINTER(c_int), POINTER(c_double), c_void_p)
peakfinder_type = ctypes.CFUNCTYPE(c_void_p, POINTER(c_int), POINTER(BOUNDS), 
  POINTER(c_int), POINTER(c_double))
  
def Vegas(integrand, ndim, userdata=NULL, 
    epsrel=EPSREL, epsabs=EPSABS, verbose=0, ncomp=1, seed=None,
    mineval=MINEVAL, maxeval=MAXEVAL, nstart=NSTART, 
    nincrease=NINCREASE, nbatch=NBATCH,
    gridno=GRIDNO, statefile=NULL):
  
  neval = c_int()
  fail = c_int()
  comp = c_int()
  ARR = c_double * NCOMP
  integral = ARR()
  error = ARR()
  prob  = ARR()
  
  if seed is None:
    seed = 0
  
  lib.Vegas(ndim, ncomp, integrand_type(integrand), userdata,
    c_double(epsrel), c_double(epsabs), verbose, seed,
    mineval, maxeval, nstart, nincrease, nbatch,
    gridno, statefile,
    byref(neval), byref(fail), integral, error, prob)
  
  return dict(neval=neval.value, fail=fail.value, comp=comp.value,
    results=[{
      'integral':integral[comp], 
      'error':error[comp], 
      'prob':prob[comp]
      } for comp in range(ncomp)])

def Suave(integrand, ndim, nnew=1000, flatness=25., userdata=NULL, 
    epsrel=EPSREL, epsabs=EPSABS, verbose=0, ncomp=1, seed=None,
    mineval=MINEVAL, maxeval=MAXEVAL):
  
  neval = c_int()
  fail = c_int()
  comp = c_int()
  nregions = c_int()
  ARR = c_double * NCOMP
  integral = ARR()
  error = ARR()
  prob  = ARR()
    
  if seed is None:
    seed = 0
  
  lib.Suave(ndim, ncomp, integrand_type(integrand), userdata,
    c_double(epsrel), c_double(epsabs), verbose, seed,
    mineval, maxeval, nnew, c_double(flatness),
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
    epsrel=EPSREL, epsabs=EPSABS, verbose=0):
  
  neval = c_int()
  fail = c_int()
  comp = c_int()
  nregions = c_int()
  ARR = c_double * NCOMP
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

  lib.Divonne(ndim, ncomp, integrand_type(integrand), userdata,
    c_double(epsrel), c_double(epsabs), verbose, seed,
    mineval, maxeval, key1, key2, key3, maxpass, 
    c_double(border), c_double(maxchisq), c_double(mindeviation), 
    ngiven, ldxgiven, xgiven, nextra, peakfinder, 
    byref(nregions), byref(neval), byref(fail), integral, error, prob)
  
  return dict(neval=neval.value, fail=fail.value, comp=comp.value, nregions=nregions.value,
    results=[{
      'integral':integral[comp], 
      'error':error[comp], 
      'prob':prob[comp]
      } for comp in range(ncomp)])

def Cuhre(integrand, ndim, 
    key, mineval=MINEVAL, maxeval=MAXEVAL, ncomp=1,
    userdata=NULL, seed=None,
    epsrel=EPSREL, epsabs=EPSABS, verbose=0):
  
  neval = c_int()
  fail = c_int()
  comp = c_int()
  nregions = c_int()
  ARR = c_double * NCOMP
  integral = ARR()
  error = ARR()
  prob  = ARR()
  
  if seed is None:
    seed = 0

  lib.Cuhre(ndim, ncomp, integrand_type(integrand), userdata,
    c_double(epsrel), c_double(epsabs), verbose,
    mineval, maxeval, key, 
    byref(nregions), byref(neval), byref(fail), integral, error, prob)
  
  return dict(neval=neval.value, fail=fail.value, comp=comp.value, nregions=nregions.value,
    results=[{
      'integral':integral[comp], 
      'error':error[comp], 
      'prob':prob[comp]
      } for comp in range(ncomp)])

if __name__ == '__main__':
  sq = lambda x: x**2
  import math

  def Integrand(ndim, xx, ncomp, ff, userdata):
    x,y,z = [xx[i] for i in range(ndim.contents.value)]
    result = math.sin(x)*math.cos(y)*math.exp(z)
    ff[0] = result
    return 0
    

  NDIM = 3
  NCOMP = 1

  NNEW = 1000
  FLATNESS = 25.

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

  KEY = 0
  
  def print_header(name):
    print '-------------------- %s test -------------------' % name
  def print_results(name, results):
    keys = ['nregions', 'neval', 'fail']
    keys = filter(results.has_key, keys)
    text = map(lambda k: "%s %d" % (k, results[k]), keys)
    print "%s RESULT:\t" % name.upper() + "\t".join(text)
    for comp in results['results']:
      print "%s RESULT:\t" % name.upper() + \
	"%(integral).8f +- %(error).8f\tp = %(prob).3f\n" % comp
  
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
  




