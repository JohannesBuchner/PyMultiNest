"""
PyAPEMoST is a module for accessing the MCMC sampler APEMoST 
(http://apemost.sf.net/)

APEMoST (Automated Parameter Estimation and Model Selection Toolkit)
is a free, fast MCMC engine that allows the user to apply 
Bayesian inference for parameter estimation and model selection.

Nota Bene: The invocation of APEMoST is not side-effect-free: 
You need to provide the usual "data" and "params" file
in the launch directory. See the APEMoST manual http://apemost.sf.net/

To change the behavior of the MCMC algorithm, recompile the APEMoST library with the 
relevant flags  (see the `API documentation <http://apemost.sourceforge.net/doc/api/html/>`_).

"""

from __future__ import absolute_import
from .run import set_function, MCMC, calibrate_first_chain, calibrate_other_chains, calibrate, run
from .analyse import create_histogram, create_histograms, print_model_probability, model_probability
from .watch import ProgressWatcher


