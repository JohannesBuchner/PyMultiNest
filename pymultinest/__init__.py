"""
PyMultiNest is a module to use the MultiNest sampling engine.

To import this module, you must have *libmultinest.so* (MultiNest) compiled 
and in your LD_LIBRARY_PATH.

  Otherwise you will get an error like this::

    > OSError: libmultinest.so: cannot open shared object file: No such file or directory

Common parameters:

outputfiles_basename is the prefix used for the output files of MultiNest 
(default chains/1-).

"""
from __future__ import absolute_import, unicode_literals, print_function
try:
	from .run import run
	from .solve import solve, Solver
except ImportError as e:
	print('WARNING: Only MultiNest analysing capabilities enabled, no running.')
	print('WARNING: check installed packages (import failed with: "%s")' % e)
from .analyse import Analyzer
try:
	from .watch import ProgressWatcher, ProgressPrinter, ProgressPlotter
	from .plot import PlotMarginal, PlotMarginalModes
except ImportError as e:
	print('WARNING: Only MultiNest run capabilities enabled, no plotting available')
	print('WARNING: check matplotlib installation (import failed with "%s")' % e)

