#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "sample_beta.h"


//main Gibbs sampler program.
  
SEXP C_GibbsFW(SEXP R_y, SEXP R_IDL, SEXP R_IDE, SEXP R_g, SEXP R_b, SEXP R_h, SEXP R_nIter, SEXP R_burnIn, SEXP R_thin, SEXP R_saveFile, SEXP R_S, SEXP R_Sg, SEXP R_Sb, SEXP R_Sh, SEXP R_df, SEXP R_dfg, SEXP R_dfb, SEXP R_dfh,SEXP R_var_e, SEXP R_var_g, SEXP R_var_b, SEXP R_var_h,SEXP R_mu,SEXP R_L, SEXP R_Linv){
//data and initial values
    PROTECT(R_y=AS_NUMERIC(R_y));
    PROTECT(R_IDL=AS_INTEGER(R_IDL));
    PROTECT(R_IDE=AS_INTEGER(R_IDE));
    PROTECT(R_g=AS_NUMERIC(R_g));
    PROTECT(R_b=AS_NUMERIC(R_b));
    PROTECT(R_h=AS_NUMERIC(R_h));
    PROTECT(R_L=AS_NUMERIC(R_L));
    PROTECT(R_Linv=AS_NUMERIC(R_Linv));
 
    
    int i,j,k; 
    int n= length(R_y);
    int ng=length(R_g);
    int nh=length(R_h);  
    double *g=(double *) R_alloc(ng*n,sizeof(double)); 
    double *b=(double *) R_alloc(ng*n,sizeof(double));
    double *h=(double *) R_alloc(nh*n,sizeof(double)); 
    for(i=0;i<ng;i++){
   		g[i]= NUMERIC_POINTER(R_g)[i];
        b[i]= NUMERIC_POINTER(R_b)[i];
    }
    for(i=0;i<nh;i++)h[i]=NUMERIC_POINTER(R_h)[i];
   
    
       
    double *y=NUMERIC_POINTER(R_y);
    int *IDL=INTEGER_POINTER(R_IDL);
    int *IDE=INTEGER_POINTER(R_IDE);
    double *L=NUMERIC_POINTER(R_L);
    double *Linv=NUMERIC_POINTER(R_Linv);
    double mu[1];mu[0]=NUMERIC_VALUE(R_mu);
    double var_e=NUMERIC_VALUE(R_var_e);
    double var_b=NUMERIC_VALUE(R_var_b);
    double var_g=NUMERIC_VALUE(R_var_g);
    double var_h=NUMERIC_VALUE(R_var_h);   
    int nIter=INTEGER_VALUE(R_nIter);
    int burnIn=INTEGER_VALUE(R_burnIn);
    int thin=INTEGER_VALUE(R_thin);
    int nSamples=nIter-burnIn+1;
    int C_IDL[n];
    int C_IDE[n];
    for(i=0;i<n;i++){C_IDL[i]=IDL[i]-1; C_IDE[i]=IDE[i]-1;}
    
    
   // int nParameters=length(R_parameters);//will include this later so we can define parameter names to print out.
    double S=NUMERIC_VALUE(R_S);
    double Sg=NUMERIC_VALUE(R_Sg);
    double Sb=NUMERIC_VALUE(R_Sb);
    double Sh=NUMERIC_VALUE(R_Sh);
    double df=NUMERIC_VALUE(R_df);
    double dfg=NUMERIC_VALUE(R_dfg);
    double dfb=NUMERIC_VALUE(R_dfb);
    double dfh=NUMERIC_VALUE(R_dfh);  
    double SS,DF;
    
    //saveFile 
    FILE *fsaveFile = fopen(CHAR(STRING_ELT(R_saveFile,0)),"w");
 	if (fsaveFile == NULL) error("Can't open input file !\n");

    

/************************************************
* posteria storage
************************************************/ 
	SEXP R_post_mu,R_post_var_g,R_post_var_b,R_post_var_h,R_post_var_e,R_post_g,R_post_b,R_post_h;
	PROTECT(R_post_mu=allocVector(REALSXP,1));
	PROTECT(R_post_var_g=allocVector(REALSXP,1));
	PROTECT(R_post_var_b=allocVector(REALSXP,1));
	PROTECT(R_post_var_h=allocVector(REALSXP,1));
	PROTECT(R_post_var_e=allocVector(REALSXP,1));
	PROTECT(R_post_g=allocVector(REALSXP,ng));
	PROTECT(R_post_b=allocVector(REALSXP,ng));
	PROTECT(R_post_h=allocVector(REALSXP,nh)); 
	double post_mu=0, post_var_e=0,post_var_g=0,post_var_b=0,post_var_h=0;
	double *post_g=NUMERIC_POINTER(R_post_g); 
	double *post_b=NUMERIC_POINTER(R_post_b); 
	double *post_h=NUMERIC_POINTER(R_post_h); 
	for(i=0;i<ng;i++){post_g[i]=0;post_b[i]=0;}
	for(i=0;i<nh;i++){post_h[i]=0;}
/************************************************
* Gibbs sampler 
************************************************/ 
	GetRNGstate();
	double *e=(double *) R_alloc(n,sizeof(double));
	double *X=(double *) R_alloc(n,sizeof(double));

	double *XL, *Xvec,*delta_g,*delta_b,*post_delta_g;

	if(!ISNAN(L[0])){
		XL=(double *) R_alloc(n*ng,sizeof(double));
		Xvec=(double *) R_alloc(n*ng,sizeof(double));
		delta_g=(double *)R_alloc(ng,sizeof(double));
		delta_b=(double *)R_alloc(ng,sizeof(double));
		post_delta_g=(double *)R_alloc(ng,sizeof(double));
		for(j=0;j<ng;j++){
			post_delta_g[j]=0;
			for(i=0;i<n;i++){
			XL[n*j+i]=L[j*ng+C_IDL[i]];
			}
		}
		Ldelta(delta_g,Linv,g,ng);
		Ldelta(delta_b,Linv,b,ng);
	}


//*initial values for e. 
	for(j=0;j<n;j++) e[j]=y[j]-mu[0]-g[C_IDL[j]]-(1+b[C_IDL[j]])*(h[C_IDE[j]]);

//begin Gibbs sampler
//no genotypic correlations
	if(ISNAN(L[0])){
		//Rprintf("start Gibbs sampler, mu:%.2f, var_e:%.2f, var_g:%.2f:, var_h:%.2f,var_b:%.2f,g[0]:%.2f,b[0]:%.2f,h[0]:%.2f\n",mu[0],var_e,var_g,var_b,var_h,g[0],b[0],h[0]);
		for(i=0; i<nIter;i++){
			//sample environment effect h
			for(j=0;j<n;j++)X[j]=(1.0+b[C_IDL[j]]);
			sample_beta_ID(h,e,C_IDE,X,n,nh,var_e,var_h);
	   
	    	//sample b and g separately;
				//sample  b
			for(j=0;j<n;j++) X[j]=h[C_IDE[j]];
			sample_beta_ID(b,e,C_IDL,X,n,ng,var_e,var_b);
				//sample  g
			sample_beta_ID_x1(g,e,C_IDL,n,ng,var_e,var_g);
		
		/*	//sample b and g together
			double tXX;
			double tXy;
			for(j=0;j<n;j++) X[j]=h[C_IDE[j]];
			for(k=0;k<ng;k++){
				tXy=0;
				tXX=0;
				//sample b
				for(j=0;j<n;j++){
					if(C_IDL[j]==k){
						e[j]=e[j]+b[C_IDL[j]]*X[j];
						tXX+=X[j]*X[j];
						tXy+=X[j]*e[j];
					}
				}
				b[k]=sample_betaj(tXX,tXy,var_e,var_b);
				for(j=0;j<n;j++){
					if(C_IDL[j]==k){
						e[j]=e[j]-b[C_IDL[j]]*X[j];
					}
				}
				// sample g
				tXy=0;
				tXX=0;
				for(j=0;j<n;j++){
					if(C_IDL[j]==k){
						e[j]=e[j]+g[C_IDL[j]];
						tXX+=1;
						tXy+=e[j];
					}
				}	
				g[k]=sample_betaj(tXX,tXy,var_e,var_g);	
				for(j=0;j<n;j++){
					if(C_IDL[j]==k){
						e[j]=e[j]-g[C_IDL[j]];
					}
				}
			}
		*/	
			//var_e;
			SS=S;
			for(j=0;j<n;j++) SS+=e[j]*e[j];
			DF=n+df;
			var_e=SS/rchisq(DF);
			//var_h
			SS=Sh;
			for(j=0;j<nh;j++) SS+=h[j]*h[j];
			DF=nh+dfh;
			var_h=SS/rchisq(DF);
			//var_b;
			SS=Sb;
			for(j=0;j<ng;j++) SS+=b[j]*b[j];
			DF=ng+dfb;
			var_b=SS/rchisq(DF);
			//var_g;
    		SS=Sg;
			for(j=0;j<ng;j++)SS+=g[j]*g[j];
			DF=ng+dfg;
			var_g=SS/rchisq(DF);

			//sample intercept
			sample_mu(mu,e,var_e,n);
        	//posteria means
			if(i>=(burnIn-1)){
				post_mu += mu[0]/nSamples;
				post_var_g += var_g/nSamples;
				post_var_h += var_h/nSamples;
				post_var_e += var_e/nSamples;
				post_var_b += var_b/nSamples;
 				for(j=0;j<ng;j++){
    				post_g[j] += g[j]/nSamples;
					post_b[j] += b[j]/nSamples;
				}
 				for(j=0;j<nh;j++) post_h[j] += h[j]/nSamples;	
			}
			if(i==0)fprintf(fsaveFile,"%s,%s,%s,%s,%s,%s,%s,%s\n","mu","var_g","var_b","var_h","var_e","g[0]","b[0]","p[0]");
			if((i+1)%thin==0)fprintf(fsaveFile,"%f,%f,%f,%f,%f,%f,%f,%f\n",mu[0],var_g,var_b,var_h,var_e,g[0],b[0],h[0]);
    		//if((i+1)%100==0){Rprintf("iter:%d\n",i+1);}//Rprintf("iter:%d, mu:%.2f, var_e:%.2f, var_g:%.2f:, var_h:%.2f,var_b:%.2f, g[0]%.2f, b[0]%.2f, h[0]%.2f\n",i+1,mu[0],var_e,var_g,var_b,var_h,g[0],b[0],h[0]);
		} //end of iterations
	}else{
	//Rprintf("start Gibbs sampler, mu:%.2f, var_e:%.2f, var_g:%.2f:, var_h:%.2f,var_b:%.2f,delta_g[0]:%.2f,b[0]:%.2f,h[0]:%.2f\n",mu[0],var_e,var_g,var_b,var_h,delta_g[0],b[0],h[0]);
		for(i=0; i<nIter;i++){
       
	  	  	//sample environment effect h
			for(j=0;j<n;j++)X[j]=(1.0+b[C_IDL[j]]);
	 		sample_beta_ID(h,e,C_IDE,X,n,nh,var_e,var_h);
	    	//sample b 
			for(j=0;j<n;j++) {
				for(k=0;k<ng;k++){
				Xvec[k*n+j]=XL[k*n+j]*h[C_IDE[j]];
				}
			}
			sample_betaX(delta_b,e,Xvec,n,ng,var_e,var_b);
			Ldelta(b,L,delta_b,ng);	
		
			//sample main effect of genotypes g
			sample_betaX(delta_g,e, XL, n, ng,var_e,var_g);	
			//var_e;
			SS=S;
			for(j=0;j<n;j++) SS+=e[j]*e[j];
			DF=n+df;
			var_e=SS/rchisq(DF);
			//var_h
			SS=Sh;
			for(j=0;j<nh;j++) SS+=h[j]*h[j];
			DF=nh+dfh;
			var_h=SS/rchisq(DF);
			//var_b;
			SS=Sb;
			for(j=0;j<ng;j++) SS+=delta_b[j]*delta_b[j];
			DF=ng+dfb;
			var_b=SS/rchisq(DF);
			//var_g;
    		SS=Sg;
			for(j=0;j<ng;j++)SS+=delta_g[j]*delta_g[j];
			DF=ng+dfg;
			var_g=SS/rchisq(DF);

		
			//sample intercept
	    	sample_mu(mu,e,var_e,n);
	    
			if(i>=(burnIn-1)){
				post_mu += mu[0]/nSamples;
				post_var_g += var_g/nSamples;
				post_var_h += var_h/nSamples;
				post_var_e += var_e/nSamples;
				post_var_b += var_b/nSamples;
 				for(j=0;j<ng;j++){
   					post_delta_g[j] += delta_g[j]/nSamples;
					post_b[j] += b[j]/nSamples;
				}
 				for(j=0;j<nh;j++) post_h[j] += h[j]/nSamples;	
			}
			if(i==0)fprintf(fsaveFile,"%s,%s,%s,%s,%s,%s,%s,%s\n","mu","var_g","var_b","var_h","var_e","delta_g[0]","b[0]","p[0]");
			if((i+1)%thin==0){    fprintf(fsaveFile,"%f,%f,%f,%f,%f,%f,%f,%f\n",mu[0],var_g,var_b,var_h,var_e,delta_g[0],b[0],h[0]);
    }
			//if((i+1)%100==0){Rprintf("iter:%d\n",i+1);}//Rprintf("iter:%d, mu:%.2f, var_e:%.2f, var_g:%.2f:, var_h:%.2f,var_b:%.2f, delta_g[0]%.2f, b[0]%.2f, h[0]%.2f\n",i+1,mu[0],var_e,var_g,var_b,var_h,delta_g[0],b[0],h[0]);}
		}//end of iterations
		Ldelta(post_g,L,post_delta_g,ng);
	}


	fclose(fsaveFile);
	PutRNGstate();//This write random numbers into R environments. 
	//store posteria means into SEXP and return

	REAL(R_post_mu)[0]=post_mu;
	REAL(R_post_var_g)[0]=post_var_g;
	REAL(R_post_var_b)[0]=post_var_b;
	REAL(R_post_var_h)[0]=post_var_h;
	REAL(R_post_var_e)[0]=post_var_e;
	SEXP list;
	PROTECT(list = allocVector(VECSXP, 8));
    SET_VECTOR_ELT(list, 0, R_post_mu);
    SET_VECTOR_ELT(list, 1, R_post_var_g);
    SET_VECTOR_ELT(list, 2, R_post_var_b);
    SET_VECTOR_ELT(list, 3, R_post_var_h);
    SET_VECTOR_ELT(list, 4, R_post_var_e);
    SET_VECTOR_ELT(list, 5, R_post_g);
    SET_VECTOR_ELT(list, 6, R_post_b);
    SET_VECTOR_ELT(list,7, R_post_h);
      
	UNPROTECT(17);

	return(list);

}	



