PyMultiNest
===========================

.. toctree::
   pymultinest_run
   pymultinest_analyse
   pymultinest_watch

.. automodule:: pymultinest
	:members:

Make sure you `installed everything <install>`_ correctly.

Demo program for PyMultiNest
--------------------------------

Take a look at the demos:

- `Demo with simple functions pymultinest_demo.py <https://github.com/JohannesBuchner/PyMultiNest/blob/master/pymultinest_demo.py>`_.
- `Demo with a class pymultinest_solver_demo.py <https://github.com/JohannesBuchner/PyMultiNest/blob/master/pymultinest_solver_demo.py>`_.
- `The old API, which allows C/Fortran Likelihood functions pymultinest_demo.py <https://github.com/JohannesBuchner/PyMultiNest/blob/master/pymultinest_demo_old.py>`_.

Tutorial
------------
A `guided tutorial <http://johannesbuchner.github.io/pymultinest-tutorial/>`_
was created for a workshop on Monte Carlo methods.

Citing PyMultiNest
-------------------------------

Please cite `MultiNest <http://ccpforge.cse.rl.ac.uk/gf/project/multinest/>`_, the 
algorithm used for computation.

If you find PyMultiNest enables your research, please consider citing my 
publication to give back for the time I invested:

`Buchner et al. 2014, A&A <http://www.aanda.org/articles/aa/abs/2014/04/aa22971-13/aa22971-13.html>`_ (BibTex there)

In this paper, I introduce the software package officially, and apply the methodology of applying MultiNest
to `X-ray spectral analysis <https://github.com/JohannesBuchner/BXA>`_.

Using PyMultiNest with MPI
----------------------------

Install the `mpi4py <http://mpi4py.scipy.org/>`_ python library. 
You can obtain it via your package manager
or via the python installer (pip or the older easy_install).

Once installed, run your python program with::

	mpiexec -n 4 python myprogram.py

PyMultiNest detects that it is run with multiple cores and will load
the "libmultinest_mpi" library instead. No further modifications are necessary.



