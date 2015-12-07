#ifndef SAMPLE_BETA_H
#define SAMPLE_BETA_H
double sample_betaj(double tXX, double tXy,double var_e,double var_b);
void sample_mu(double *mu,double *e, double var_e, int n);
void sample_beta_ID_x1(double *b, double *e, const int *C_ID, int n, int ngroups,double var_e,double var_b);
void sample_beta_ID(double *b,double *e, const int *C_ID, const double *x,  int n, int ngroups,double var_e,double var_b);
void sample_betaX(double *b,double *e, double *Xvec, int nrow, int ncol,double var_e,double var_b);
void Ldelta(double *b,const double *L, const double *delta, const int ngroups);
void Udelta(double *b, const double *U, const double *delta, const int nb, const int ndelta);
#endif

double invNormpdf(double d, double mu, double sigma,int lower_tail);
void sample_mu_slicer(double *mu, double lmu, double rmu, double *e, double var_e, int n);

