"""
PyMultiNest is a module to use the MultiNest sampling engine.

To import this module, you must have 

1. *libcnest.so* (multinest_bridge) compiled and in your LD_LIBRARY_PATH

  Otherwise you will get an error like this::

    > OSError: libcnest.so: cannot open shared object file: No such file or directory

2. *libnest3.so* (MultiNest) compiled and in your LD_LIBRARY_PATH

  Otherwise you will get an error like this::

    > OSError: libnest3.so: cannot open shared object file: No such file or directory

Common parameters:

outputfiles_basename is the prefix used for the output files of MultiNest 
(default chains/1-).

"""
from run import run
from analyse import Analyzer
try:            
	from watch import ProgressWatcher, ProgressPrinter, ProgressPlotter
	from plot import PlotMarginal, PlotMarginalModes
except ImportError as e:
	print e
	print 'WARNING: no plotting available -- check matplotlib installation and error above'
	print 'Only MultiNest run capabilities enabled.'

