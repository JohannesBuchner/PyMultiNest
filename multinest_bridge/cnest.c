#include<stdlib.h>
#include<string.h>

extern void __nested_MOD_nestrun(
	int *mmodal, int *ceff, int *nlive, double *tol, double *efr, int *ndims,
	int *nPar, int *nClsPar, int *maxModes, int *updInt, double *Ztol, 
	char *root, int *seed, int *pWrap, int *fb, int *resume, 
	void (*Loglike)(double *Cube, int *n_dim, int *n_par, double *lnew), 
	int *context);

struct problem {
	void (*Prior)(double *Cube, int n_dim, int n_par);
	double (*LogLike)(double *Cube, int n_dim, int n_par);
} p;

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

void run(
	int mmodal, int ceff, int nlive, double tol, double efr, int ndims,
	int nPar, int nClsPar, int maxModes, int updInt, double Ztol, 
	char *rootstr, int seed, int * pWrap, int fb, int resume, int context
)
{
	char root[100];
	/* filling root with spaces, because fortran likes that */
	strcpy(root, rootstr);
	memset(root + strlen(root), ' ', 100 - strlen(root));

	/* running MultiNest */
	__nested_MOD_nestrun(&mmodal, &ceff, &nlive, &tol, &efr, &ndims, 
		&nPar, &nClsPar, &maxModes, &updInt, &Ztol,
		root, &seed, pWrap, &fb, &resume, _LogLike, &context);
}

