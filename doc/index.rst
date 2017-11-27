.. pymultinest documentation master file, created by
   sphinx-quickstart on Thu Jan  5 04:13:31 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pymultinest's documentation!
=======================================

About
------------------------------

This is the documentation for two python modules for Bayesian analysis,
parameter estimation and model selection: pymultinest and pycuba.

* :doc:`PyMultiNest <pymultinest>` interacts with **MultiNest**, a Nested Sampling Monte Carlo library.
* :doc:`PyCuba <pycuba>` interacts with **cuba**, a library for various numerical integration methods.

Get PyMultiNest
-------------------------------

* Follow the instructions in the :doc:`installation guide <install>`.

* You can download the source code from the `GitHub code repository <https://github.com/JohannesBuchner/PyMultiNest>`_.

Citing PyMultiNest
-------------------------------

Please cite `MultiNest <pymultinest>` and `Cuba <pycuba>` accordingly, 
depending on which algorithm you connect your Python program to using this package.

If you find PyMultiNest enables your research, please consider citing my 
publication to give back for the time I invested:

`Buchner et al. 2014, A&A <http://www.aanda.org/articles/aa/abs/2014/04/aa22971-13/aa22971-13.html>`_ (`get BibTex <http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2014A%26A...564A.125B&data_type=BIBTEX&db_key=AST&nocookieset=1>`_)

In this paper, I introduce the software package officially, and apply the methodology of applying MultiNest
to `X-ray spectral analysis <https://github.com/JohannesBuchner/BXA>`_.


Documentation:
-------------------------------

.. toctree::
   install
   pymultinest
   pycuba
   :maxdepth: -1


Indices and tables
-------------------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

