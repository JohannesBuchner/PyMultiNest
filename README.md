PyMultiNest -- Python interface for MultiNest -- Michele's fork
===============================================================

**Note: not completely committed yet... do not use for the moment!**

[PyMultiNest](http://johannesbuchner.github.com/PyMultiNest) provides programmatic access to [MultiNest](http://ccpforge.cse.rl.ac.uk/gf/project/multinest) from Python. It does so by building a Python extension that contains a C wrapper (the "bridge") to MultiNest's Fortran. By contrast, my approach is write a C-compatible wrapper in Fortran (using the [iso_c_binding](http://fortran90.org/src/best-practices.html) interface, a Fortran 2003 feature supported by modern compilers), and include it among the MultiNest source files. PyMultiNest can then access MultiNest directly using `ctypes`.

The C interface
---------------

The single MultiNest procedure that is exposed to C (for MultiNest 3.2) is

    void run(bool nest_IS,bool nest_mmodal,bool nest_ceff, \
             int nest_nlive,double nest_tol,double nest_ef,int nest_ndims,int nest_totPar,int nest_nCdims,int maxClst, \
             int nest_updInt,double nest_Ztol,char nest_root[],int seed,int nest_pWrap[], \
             bool nest_fb,bool nest_resume,bool nest_outfile,bool initMPI,double nest_logZero,int nest_maxIter, \
             double (*loglike)(double *,int,int,void *context), \
             void (*dumper)(int,int,int,double *,double *,double *,double,double,double,void *context),void *context);

(compare with the Fortran call as documented in the MultiNest `README`) and the two callbacks are

    double loglike(double *Cube,int n_dim,int nPar,void *context);

    void dumper(int nSamples,int nlive,int nPar, \
                double *physLive,double *posterior,double *paramConstr, \
                double maxLogLike,double logZ,double logZerr,void *context);
