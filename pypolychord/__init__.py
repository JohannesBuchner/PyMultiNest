"""
PyPolyChord is a module to use the PolyChord sampling engine.

To import this module, you must have *libchord.so* (PolyChord) compiled 
and in your LD_LIBRARY_PATH.

  Otherwise you will get an error like this::

    > OSError: libchord.so: cannot open shared object file: No such file or directory

Common parameters:

outputfiles_basename is the prefix used for the output files of PolyChord 
(default chains/1-).

"""
from __future__ import absolute_import, unicode_literals, print_function
try:
	from .run import run
except ImportError as e:
	print(e)
	print('WARNING: no running available -- check library and error above')
	print('Only MultiNest analysing capabilities enabled.')
