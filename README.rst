PyMultiNest -- Python interface for MultiNest
==============================================

This library provides programmatic access to MultiNest.

What is MultiNest?
-------------------

MultiNest is a program and a sampling technique. As a Bayesian inference technique,
it allows parameter estimation and model selection. (find out more in the 
MultiNest paper, http://arxiv.org/abs/0809.3437, or in a classic MCMC sampler, 
http://apemost.sf.net/ ).

The efficient Monte Carlo algorithm for sampling the parameter space is based 
on nested sampling and the idea of disjoint multi-dimensional ellipse sampling.

MultiNest is a great and fast algorithm, but I just wish it were easier to 
access. Fortran is too difficult to read, and rewriting the code in a modern 
language would be good -- for scientific verification alone --
but is just too difficult for me. 

For the scientific community, where Python is the new lingua franca (luckily),
I provide an interface to the great MultiNest library.

What does PyMultiNest do?
--------------------------

  * Provides an easy-to-use interface to MultiNest

  * Provides integration with your existing scientific python code (numpy, scipy)

  * Allows you to write Prior & LogLikelihood functions in Python. This is 
    usally fast enough -- if not, see below.

  * (Planned) plotting and visualization of progress.

  * (Planned) easy plotting and visualization of results. (Code welcome!)

Code contributions are welcome! Contact me (buchner.johannes [ät] gmx.at).

How can I use MultiNest with Python?
--------------------------------------------
Look at the documentation available at http://johannesbuchner.github.com/PyMultiNest/doc/

What is PyAPEMoST?
--------------------------------------------
Similarly to PyMultiNest, it is an access module for a Bayesian inference engine.
However, APEMoST is a Markov Chain Monte Carlo engine. See the `documentation <http://johannesbuchner.github.com/PyMultiNest/doc/html/pyapemost>`_.

Q: Python callback functions are too slow!
-------------------------------------------
If you really identified that your callback functions are too slow, even
when using the usual tricks (numpy, etc.), you can just program them into
cnest.c, effectively making them part of the cnest library.

You still have the neat python interface (default parameters, etc.), but
achieve full execution speed, as only native code is executed while
MultiNest runs.



