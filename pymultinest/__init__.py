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
except ImportError as e:
	print(e)
	print('WARNING: no running available -- check library and error above')
	print('Only MultiNest analysing capabilities enabled.')
from .analyse import Analyzer
try:
	from .watch import ProgressWatcher, ProgressPrinter, ProgressPlotter
	from .plot import PlotMarginal, PlotMarginalModes
except ImportError as e:
	print(e)
	print('WARNING: no plotting available -- check matplotlib installation and error above')
	print('Only MultiNest run capabilities enabled.')

