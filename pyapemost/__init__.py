"""
PyAPEMoST is a module for accessing the MCMC sampler APEMoST 
(http://apemost.sf.net/)

APEMoST (Automated Parameter Estimation and Model Selection Toolkit)
is a free, fast MCMC engine that allows the user to apply 
Bayesian inference for parameter estimation and model selection.



"""


from run import set_function, MCMC, calibrate_first_chain, calibrate_other_chains, calibrate, run
from analyse import create_histogram, create_histograms, print_model_probability, model_probability
from watch import ProgressWatcher


