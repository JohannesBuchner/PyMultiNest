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


Q: Python callback functions are too slow!
-------------------------------------------
If you really identified that your callback functions are too slow, even
when using the usual tricks (numpy, etc.), you can just program them into
cnest.c, effectively making them part of the cnest library.

You still have the neat python interface (default parameters, etc.), but
achieve full execution speed, as only native code is executed while
MultiNest runs.







Q: libcnest.so doesn't load
-----------------------------------------------------

Make sure it exists. Tell the program where it is by including the current 
directory into the LD_LIBRARY_PATH.

   export LD_LIBRARY_PATH=/usr/lib:/lib:.


If you receive the error that libnest3.so couldn't be loaded, 
make sure that the libnest3.so you created from MultiNest is there, and it 
really is a shared library.

$ ldd ../libnest3.so 
	linux-gate.so.1 =>  (0x0095b000)
	liblapack.so.3 => /usr/lib/atlas/liblapack.so.3 (0x00110000)
	libpthread.so.0 => /lib/libpthread.so.0 (0x006a5000)
	libgfortran.so.3 => /usr/lib/libgfortran.so.3 (0x00b63000)
	libm.so.6 => /lib/libm.so.6 (0x006bf000)
	libgcc_s.so.1 => /lib/libgcc_s.so.1 (0x006e9000)
	libquadmath.so.0 => /usr/lib/libquadmath.so.0 (0x00e70000)
	libc.so.6 => /lib/libc.so.6 (0x0095c000)
	libf77blas.so.3 => /usr/lib/atlas/libf77blas.so.3 (0x00e2a000)
	libcblas.so.3 => /usr/lib/atlas/libcblas.so.3 (0x00706000)
	/lib/ld-linux.so.2 (0x42d79000)
	libatlas.so.3 => /usr/lib/atlas/libatlas.so.3 (0x00ee3000)


This can be done, for instance, by running

gfortran -shared -llapack -lpthread -o  libnest3.so utils.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o 
