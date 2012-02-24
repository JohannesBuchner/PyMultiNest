#include<stdlib.h>
#include<stdio.h>
#include<string.h>


#define LOGLIKETYPE(f) double (f)(double *Cube, int n_dim, int n_par)
#define PRIORTYPE(f) void (f)(double *Cube, int n_dim, int n_par)

#define MULTINEST_CALLBACK(f) void (f) (double *Cube, int *ndim, int *npars, double *lnew)
#define MULTINEST_DUMPERTYPE(f) void (f)(int *nsamples, int *nlive, int *npar, \
		double **physlive, double **posterior, double **paramconstr, \
		double * maxloglike, double * logz, double *logzerr)
#define DUMPERTYPE(f) void (f)(int nsamples, int nlive, int npar, \
		double ** physlive, double ** posterior, \
		double *mean, double *std, double *best, double *map, \
		double maxloglike, double logz, double logzerr)


extern void MULTINEST_CALL(
	int *mmodal, int *ceff, int *nlive, double *tol, double *efr, int *ndims,
	int *nPar, int *nClsPar, int *maxModes, int *updInt, double *Ztol, 
	char *root, int *seed, int *pWrap, int *fb, int *resume,
	int *outfile, int *initMPI, double *nestlogzero,  
	MULTINEST_CALLBACK(*Loglike), 
	MULTINEST_DUMPERTYPE(*Dumper),
	int *context);

struct problem {
	PRIORTYPE(*Prior);
	LOGLIKETYPE(*LogLike);
	DUMPERTYPE(*Dumper);
} p = {NULL, NULL, NULL};

MULTINEST_CALLBACK(_LogLike)
{
	p.Prior(Cube, *ndim, *npars);
	*lnew = p.LogLike(Cube, *ndim, *npars);
}
MULTINEST_DUMPERTYPE(_noop_dumper)
{
}
MULTINEST_DUMPERTYPE(_DumperConverter)
{
        int i, j;
	printf("\tconverting for access ...\n");
        int n = *npar;

        double postdist[*nsamples][n + 2];
        double pLivePts[*nlive][n + 1];
        double *mean = &(paramconstr[0][0]);
        double *std = &(paramconstr[0][n]);
        double *best = &(paramconstr[0][n*2]);
        double *map = &(paramconstr[0][n*3]);
        
        printf("\t--> DumperConverter called with: %d %d %d\n", *nsamples, *nlive, *npar);
	
	for( j = 0; j < *nsamples; j++ )
		for( i = 0; i < n + 2; i++ ) // last is likelihood value & posterior
			postdist[j][i] = posterior[0][i * (*nsamples) + j];

	for( j = 0; j < *nlive; j++ )
		for( i = 0; i < n + 1; i++ ) // last is likelihood value
			pLivePts[j][i] = physlive[0][i * (*nlive) + j];

	printf("\tconverting for access ... done. \n");
	printf("\tdumper set to %p\n", p.Dumper);
	if (p.Dumper != NULL) {
		printf("\tcalling client. \n");
		printf("\tcalling %p\n", p.Dumper);
		p.Dumper(*nsamples, *nlive, n, 
			NULL /*, (double**)postdist, 
			mean, std, best, map, *maxloglike, *logz, *logzerr */);
	}
}

void reset() {
	p.Prior = NULL;
	p.LogLike = NULL;
	p.Dumper = NULL;
}

void set_function(LOGLIKETYPE(*LogLike), PRIORTYPE(*Prior)) {
	p.Prior = Prior;
	p.LogLike = LogLike;
}

void set_dumper(DUMPERTYPE(*Dumper)) {
	p.Dumper = Dumper;
}

void run(
	int mmodal, int ceff, int nlive, double tol, double efr, int ndims,
	int nPar, int nClsPar, int maxModes, int updInt, double Ztol, 
	char *rootstr, int seed, int * pWrap, int fb, int resume, int outfile, 
	int initMPI, double logZero, int context
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
	DUMPERTYPE(*dumpfunc) = _DumperConverter;
	if (p.Dumper == NULL) 
		dumpfunc = _noop_dumper;
	
	
	/* running MultiNest */
	MULTINEST_CALL(&mmodal, &ceff, &nlive, &tol, &efr, &ndims, 
		&nPar, &nClsPar, &maxModes, &updInt, &Ztol,
		root, &seed, pWrap, &fb, &resume,
		&outfile, &initMPI, &logZero,  
			 _LogLike, dumpfunc, &context);
}

