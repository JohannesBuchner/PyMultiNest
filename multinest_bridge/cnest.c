#include<stdlib.h>
#include<stdio.h>
#include<string.h>

extern void __nested__nestrun(
	int *mmodal, int *ceff, int *nlive, double *tol, double *efr, int *ndims,
	int *nPar, int *nClsPar, int *maxModes, int *updInt, double *Ztol, 
	char *root, int *seed, int *pWrap, int *fb, int *resume,
	int *outfile, int *initMPI, double *nestlogzero,  
	void (*Loglike)(double *Cube, int *n_dim, int *n_par, double *lnew), 
	void (*dumper)(int * nsamples, int *nlive, int *npar, 
		double * physlive, double *posterior, double *paramconstr,
		double * maxloglike, double * logz, double *logzerr),
	int *context);

struct problem {
	void (*Prior)(double *Cube, int n_dim, int n_par);
	double (*LogLike)(double *Cube, int n_dim, int n_par);
} p = {NULL, NULL};

void _LogLike(double *Cube, int *ndim, int *npars, double *lnew)
{
	p.Prior(Cube, *ndim, *npars);
	*lnew = p.LogLike(Cube, *ndim, *npars);
}

void set_function(
	void (*Prior)(double *Cube, int n_dim, int n_par),
        double (*LogLike)(double *Cube, int n_dim, int n_par)
) {
	p.Prior = Prior;
	p.LogLike = LogLike;
}

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *logZerr)
{
	// convert the 2D Fortran arrays to C arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;
	
	double postdist[*nSamples][*nPar + 2];
	for( i = 0; i < *nPar + 2; i++ )
		for( j = 0; j < *nSamples; j++ )
			postdist[j][i] = posterior[0][i * (*nSamples) + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[*nlive][*nPar + 1];
	for( i = 0; i < *nPar + 1; i++ )
		for( j = 0; j < *nlive; j++ )
			pLivePts[j][i] = physLive[0][i * (*nlive) + j];
}

void run(
	int mmodal, int ceff, int nlive, double tol, double efr, int ndims,
	int nPar, int nClsPar, int maxModes, int updInt, double Ztol, 
	char *rootstr, int seed, int * pWrap, int fb, int resume, int context
)
{
	int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
	double logZero = -1e100;			// points with loglike < logZero will be ignored by MultiNest
	int outfile = 1;				// write output files?

	char root[100];
	/* filling root with spaces, because fortran likes that */
	strcpy(root, rootstr);
	memset(root + strlen(root), ' ', 100 - strlen(root));
	
	if (p.Prior == NULL || p.LogLike == NULL) {
		fprintf(stderr, "Need to call set_function(prior, loglike) first, to define callback functions!\n");
		return;
	}
	
	/* running MultiNest */
	__nested__nestrun(&mmodal, &ceff, &nlive, &tol, &efr, &ndims, 
		&nPar, &nClsPar, &maxModes, &updInt, &Ztol,
		root, &seed, pWrap, &fb, &resume,
		&outfile, &initMPI, &logZero,  
			 _LogLike, dumper, &context);
}

