//use Makefile to compile this code
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "sample_beta.h"

double sample_betaj(double tXX, double tXy,double var_e,double var_b){
double C = tXX/var_e+1/var_b;
double rhs= tXy/var_e;
double b=rhs/C+sqrtf(1/C)*norm_rand();
return(b);
}

void sample_mu(double *mu,double *e, double var_e, int n){
//update mu and e.
double muhat, V_muhat;
int i;
for(i=0;i<n;i++) e[i]=e[i]+mu[0];
muhat=0;
for(i=0;i<n;i++){muhat+=e[i];}
muhat=muhat/n;
V_muhat=var_e/n;
mu[0]=muhat+sqrtf(V_muhat)*norm_rand();
for(i=0;i<n;i++) e[i]=e[i]-mu[0];
}

void sample_beta_ID_x1(double *b, double *e, const int *C_ID, int n, int ngroups,double var_e,double var_b){
//update b and e.
//x are all ones in this case y=b[ID]
//sample beta with group ID. (C_ID is the integer ID in the index convention of C). Every observation can have only one group ID
//b is the holder for the sampled beta.
//n size of ID, y
//ngroups: number of grouping factors
// var_e: residual variance
// var_b variance for the regression coefficient
int i,j;
double tXy[ngroups]; //must initialize tXy and tXX, otherwise, something unpredictable will occur
double tXX[ngroups];
for(j=0;j<n;j++) e[j]=e[j]+b[C_ID[j]];

for(i=0;i<ngroups;i++){
tXy[i]=0;
tXX[i]=0;
}

for(i=0;i<n;i++){
j=C_ID[i];
tXy[j]+=e[i];
tXX[j]+=1;
}
for(j=0;j<ngroups;j++) b[j]=sample_betaj(tXX[j],tXy[j],var_e,var_b);

for(j=0;j<n;j++) e[j]=e[j]-b[C_ID[j]];

}

//sample beta based on x, y and grouping factors(ID).
void sample_beta_ID(double *b,double *e, const int *C_ID, const double *x,  int n, int ngroups,double var_e,double var_b){
//update b and e
//e=x#b[ID]
//sample beta with group ID. (C_ID is the integer ID in the index convention of C). Every observation can have only one group ID
//b is the holder for the sampled beta.
//n size of ID, X and y
//ngroups: number of grouping factors
// var_e: residual variance
// var_b variance for the regression coefficient
int i,j;
double tXy[ngroups];
double tXX[ngroups];
for(i=0;i<ngroups;i++){
tXy[i]=0;
tXX[i]=0;
}
for(j=0;j<n;j++) e[j]=e[j]+b[C_ID[j]]*x[j];
for(i=0;i<n;i++){
j=C_ID[i];
tXy[j]+=x[i]*e[i];
tXX[j]+=x[i]*x[i];
}
for(j=0;j<ngroups;j++) b[j]=sample_betaj(tXX[j],tXy[j],var_e,var_b);
for(j=0;j<n;j++) e[j]=e[j]-b[C_ID[j]]*x[j];
}



void sample_betaX(double *b,double *e, double *Xvec, int nrow, int ncol,double var_e,double var_b){
//b and e will be updated
//sample beta from incidence matrix X. All betaj have weights in each individual, so the total non-zero values in incidence matrix X is nrow*ncol.
//b is the holder for the sampled beta.
//Xvec: incidence matrix X transformed into matrix
//nrow :numer of rows for incidence matrix C size of ID, X and y
//ncol: number of columns in incidence matrix X
// var_e: residual variance
// var_b variance for the regression coefficient
double tXy;
double tXX;
int i,j;
double *xj;

for(j=0; j<ncol;j++)
    {
	  tXy=0;
	  tXX=0;
	  xj=Xvec+j*nrow;
	    for(i=0; i<nrow; i++)
	    {
	      e[i] = e[i] + b[j]*xj[i];//e=e+bj*xj
	      tXy+=xj[i]*e[i];
	      tXX+=xj[i]*xj[i];
	    }
	  b[j]=sample_betaj(tXX, tXy,var_e,var_b);
	  
	  for(i=0; i<nrow; i++)
	  {
	    e[i] = e[i] - b[j]*xj[i];
	  }
	  }
}

//calculate b=L%*%delta
void Ldelta(double *b,const double *L, const double *delta, const int ngroups){
int i,j;
for(i=0;i<ngroups;i++) {
 b[i]=0; 
for(j=0;j<=i;j++){
 b[i]+=L[j*ngroups+i]*delta[j];
}
}
}

