#include<stdlib.h>
#include<stdio.h>
#include<string.h>


#define LOGLIKETYPE(f) double (*f)(double *Cube, int n_dim, int n_par)
#define PRIORTYPE(f) void (*f)(double *Cube, int n_dim, int n_par)

#define MULTINEST_CALLBACK(f) void f (double *Cube, int *ndim, int *npars, double *lnew)
#define DUMPERTYPE(f) void (*f)(int *nsamples, int *nlive, int *npar, \
		double **physlive, double **posterior, double **paramconstr, \
		double * maxloglike, double * logz, double *logzerr)


extern void MULTINEST_CALL(
	int *mmodal, int *ceff, int *nlive, double *tol, double *efr, int *ndims,
	int *nPar, int *nClsPar, int *maxModes, int *updInt, double *Ztol, 
	char *root, int *seed, int *pWrap, int *fb, int *resume,
	int *outfile, int *initMPI, double *nestlogzero,  
	MULTINEST_CALLBACK(Loglike), 
	DUMPERTYPE(Dumper),
	int *context);

struct problem {
	PRIORTYPE(Prior);
	LOGLIKETYPE(LogLike);
	DUMPERTYPE(Dumper);
} p = {NULL, NULL, NULL};

MULTINEST_CALLBACK(_LogLike)
{
	p.Prior(Cube, *ndim, *npars);
	*lnew = p.LogLike(Cube, *ndim, *npars);
}

void set_function(LOGLIKETYPE(LogLike), PRIORTYPE(Prior)) {
	p.Prior = Prior;
	p.LogLike = LogLike;
}

void set_dumper(DUMPERTYPE(Dumper)) {
	p.Dumper = Dumper;
}

void run(
	int mmodal, int ceff, int nlive, double tol, double efr, int ndims,
	int nPar, int nClsPar, int maxModes, int updInt, double Ztol, 
	char *rootstr, int seed, int * pWrap, int fb, int resume, int context, 
	int initMPI, double logZero, int outfile
)
{
	char root[100];
	/* filling root with spaces, because fortran likes that */
	strcpy(root, rootstr);
	memset(root + strlen(root), ' ', 100 - strlen(root));
	
	if (p.Prior == NULL || p.LogLike == NULL) {
		fprintf(stderr, "Need to call set_function(prior, loglike) first, to define callback functions!\n");
		return;
	}
	
	/* running MultiNest */
	MULTINEST_CALL(&mmodal, &ceff, &nlive, &tol, &efr, &ndims, 
		&nPar, &nClsPar, &maxModes, &updInt, &Ztol,
		root, &seed, pWrap, &fb, &resume,
		&outfile, &initMPI, &logZero,  
			 _LogLike, p.Dumper, &context);
}

