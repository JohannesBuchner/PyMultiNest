PyMultiNest -- Python interface for MultiNest
==============================================

This library provides programmatic access to MultiNest and PyCuba.

What is MultiNest?
-------------------

MultiNest is a program and a sampling technique. As a Bayesian inference technique,
it allows parameter estimation and model selection. (find out more in the 
MultiNest paper, http://arxiv.org/abs/0809.3437, or in a classic MCMC sampler, 
http://apemost.sf.net/ ). Recently, MultiNest added Importance Nested Sampling 
(INS, see http://arxiv.org/abs/1306.2144) which is now also supported.

The efficient Monte Carlo algorithm for sampling the parameter space is based 
on nested sampling and the idea of disjoint multi-dimensional ellipse sampling.

For the scientific community, where Python is becoming the new lingua franca (luckily),
I provide an interface to MultiNest.

What does PyMultiNest do?
--------------------------

PyMultiNest 

  * provides an easy-to-use interface to MultiNest

  * provides integration with your existing scientific Python code (numpy, scipy)

  * allows you to write Prior & LogLikelihood functions in Python.

PyMultiNest can 

  * Plot and visualize MultiNests progress (watch.ProgressWatcher, watch.ProgressPlotter). This is still fairly basic, contributions and ideas are welcome)

  * Easy plotting, visualization and summary of MultiNest results.

The plotting can be run on existing MultiNest output, and when not using PyMultiNest for running MultiNest.

Code contributions are welcome! Contact me (buchner.johannes [ät] gmx.at).

How can I use MultiNest with Python?
--------------------------------------------
Look at the documentation available at http://johannesbuchner.github.com/PyMultiNest/index.html

Citing PyMultiNest
--------------------------------------------
See the documentation at http://johannesbuchner.github.com/PyMultiNest/index.html

What is PyAPEMoST?
--------------------------------------------
Similarly to PyMultiNest, it is an access module for a Bayesian inference engine.
However, APEMoST is a Markov Chain Monte Carlo engine. See the `documentation <http://johannesbuchner.github.com/PyMultiNest/pyapemost>`_.

What is PyCuba?
--------------------------------------------
Cuba (http://www.feynarts.de/cuba/, https://github.com/JohannesBuchner/cuba) is a multidimensional numerical integration library for low dimensions. PyCuba allows integration of Python functions, providing an advanced alternative to the basic functions provided in scipy.integrate.

In the Bayesian sense, it is possible to use Cuba for model selection.

Q: Python callback functions are too slow!
-------------------------------------------
If you really identified that your callback functions are too slow, even
when using the usual tricks (numpy, etc.), you can implement and compile 
them as C functions.

You still have the neat python interface (default parameters, etc.), but
achieve full execution speed, as only native code is executed while
MultiNest runs.



