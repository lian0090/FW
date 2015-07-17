#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "sample_beta.h"
#include <stdlib.h>
#include <string.h>

//create a new string, which is the combination of two strings.
char* concat(const char *s1, char *s2)
{
    char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the zero-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}


//main Gibbs sampler program.

SEXP C_GibbsFW(SEXP R_y, SEXP R_IDL, SEXP R_IDE, SEXP R_g, SEXP R_b, SEXP R_h, SEXP R_nIter, SEXP R_burnIn, SEXP R_thin, SEXP R_saveFile, SEXP R_S, SEXP R_Sg, SEXP R_Sb, SEXP R_Sh, SEXP R_df, SEXP R_dfg, SEXP R_dfb, SEXP R_dfh,SEXP R_var_e, SEXP R_var_g, SEXP R_var_b, SEXP R_var_h,SEXP R_mu,SEXP R_L, SEXP R_Linv, SEXP R_whNA , SEXP R_whNotNA, SEXP R_VARstore, SEXP R_ENVstore){
    int nProtect=0;
    //data and initial values
    PROTECT(R_y=AS_NUMERIC(R_y));  nProtect+=1;
    PROTECT(R_IDL=AS_INTEGER(R_IDL)); nProtect+=1;
    PROTECT(R_IDE=AS_INTEGER(R_IDE)); nProtect+=1;
    PROTECT(R_g=AS_NUMERIC(R_g)); nProtect+=1;
    PROTECT(R_b=AS_NUMERIC(R_b)); nProtect+=1;
    PROTECT(R_h=AS_NUMERIC(R_h)); nProtect+=1;
    PROTECT(R_L=AS_NUMERIC(R_L)); nProtect+=1;
    PROTECT(R_Linv=AS_NUMERIC(R_Linv)); nProtect+=1;
    PROTECT(R_whNA=AS_INTEGER(R_whNA)); nProtect+=1;
    PROTECT(R_VARstore=AS_INTEGER(R_VARstore)); nProtect+=1;
    PROTECT(R_ENVstore=AS_INTEGER(R_ENVstore)); nProtect+=1;   
    int i,j,k;
    int n= length(R_y);
    int ng=length(R_g);
    int nh=length(R_h);
    int nNa=length(R_whNA);
    int nNotNa=length(R_whNotNA);
    int nVAR_Store=length(R_VARstore);
    int nENV_Store=length(R_ENVstore);
    double mu[1];mu[0]=NUMERIC_VALUE(R_mu);
    double var_e=NUMERIC_VALUE(R_var_e);
    double var_b=NUMERIC_VALUE(R_var_b);
    double var_g=NUMERIC_VALUE(R_var_g);
    double var_h=NUMERIC_VALUE(R_var_h);
    int nIter=INTEGER_VALUE(R_nIter);
    int burnIn=INTEGER_VALUE(R_burnIn);
    int thin=INTEGER_VALUE(R_thin);
    int nSamples=(nIter-burnIn)/thin;//burnIn is the number of discarded samples
    double *y=NUMERIC_POINTER(R_y);
    int *whNA=INTEGER_POINTER(R_whNA);
    int *whNotNA=INTEGER_POINTER(R_whNotNA);
    int *IDL=INTEGER_POINTER(R_IDL);
    int *IDE=INTEGER_POINTER(R_IDE);
    double *L=NUMERIC_POINTER(R_L);
    double *Linv=NUMERIC_POINTER(R_Linv);
    int *VARstore=INTEGER_POINTER(R_VARstore);
    int *ENVstore=INTEGER_POINTER(R_ENVstore);
    
    int C_IDL[n];
    int C_IDE[n];
    for(i=0;i<n;i++){C_IDL[i]=IDL[i]-1; C_IDE[i]=IDE[i]-1;}
  
 

    
    //g, b, h are duplicates of R_g, R_b, R_h, they do not point to R_g, R_b, R_h, and do not modify them
    double *g=(double *) R_alloc(ng,sizeof(double));
    double *b=(double *) R_alloc(ng,sizeof(double));
    double *h=(double *) R_alloc(nh,sizeof(double));
    double *yhat=(double *) R_alloc(n,sizeof(double));
    double *yStar=(double *)R_alloc(n,sizeof(double));
    for(i=0;i<ng;i++){
        g[i]= NUMERIC_POINTER(R_g)[i];
        b[i]= NUMERIC_POINTER(R_b)[i];
    }
    for(i=0;i<nh;i++)h[i]=NUMERIC_POINTER(R_h)[i];
    for(i=0;i<n;i++){
    yhat[i]=mu[0];
    yStar[i]=y[i];
    }
    //starting value for yStar[whNA] is set to mu[0]
     if(nNa>0){
            for(j=0;j<nNa;j++){
                yStar[(whNA[j]-1)]=mu[0];
            }
        }
    
    
    
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
   //The following code are to save each parameter separately in a folder.
   // const char *saveAt =CHAR(STRING_ELT(R_saveAt, 0));
   // FILE *fmu=fopen(concat(saveAt,"mu.dat"));
   // FILE *fvar_e=fopen(concat(saveAt,"var_e.dat"));
    //FILE *fvar_g=fopen(concat(saveAt,"var_g.dat"));
    //FILE *fvar_b=fopen(concat(saveAt,"var_b.dat"));
    //FILE *fvar_h=fopen(concat(saveAt,"var_h.dat"));
    //FILE *fg=fopen(concat(saveAt,"g.dat"));
    //FILE *fb=fopen(concat(saveAt,"b.dat"));
    //FILE *fh=fopen(concat(saveAt,"h.dat"));
    
    //headers for samples file.
  
    fprintf(fsaveFile,"%s,%s,%s,%s,%s","mu","var_g","var_b","var_h","var_e");
            if(ISNAN(L[0])){
            	for(j=0;j<nVAR_Store;j++){
            	fprintf(fsaveFile,",g[%d]",VARstore[j]);
            	}
            }else{
            	for(j=0;j<nVAR_Store;j++){
            	fprintf(fsaveFile,",delta_g[%d]",VARstore[j]);
            	}
            }	
            for(j=0;j<nVAR_Store;j++){
            	fprintf(fsaveFile,",b[%d]",VARstore[j]);
            }
            for(j=0;j<nENV_Store;j++){
            	fprintf(fsaveFile,",h[%d]",ENVstore[j]);
            }
            
    fprintf(fsaveFile,"\n");
            
        
    /************************************************
     * posteria and yhat storage
     ************************************************/
    SEXP R_post_mu,R_post_var_g,R_post_var_b,R_post_var_h,R_post_var_e,R_post_g,R_post_b,R_post_h, R_post_yhat, R_post_logLik;
    PROTECT(R_post_mu=allocVector(REALSXP,1)); nProtect+=1;
    PROTECT(R_post_var_g=allocVector(REALSXP,1)); nProtect+=1;
    PROTECT(R_post_var_b=allocVector(REALSXP,1));nProtect+=1;
    PROTECT(R_post_var_h=allocVector(REALSXP,1));nProtect+=1;
    PROTECT(R_post_var_e=allocVector(REALSXP,1));nProtect+=1;
    PROTECT(R_post_logLik=allocVector(REALSXP,1));nProtect+=1;
    PROTECT(R_post_g=allocVector(REALSXP,ng));nProtect+=1;
    PROTECT(R_post_b=allocVector(REALSXP,ng));nProtect+=1;
    PROTECT(R_post_h=allocVector(REALSXP,nh));nProtect+=1;
    PROTECT(R_post_yhat=allocVector(REALSXP,n));nProtect+=1;
   
   

    
    double post_mu=0, post_var_e=0,post_var_g=0,post_var_b=0,post_var_h=0,post_logLik=0,logLik=0;
    double *post_g=NUMERIC_POINTER(R_post_g);
    double *post_b=NUMERIC_POINTER(R_post_b);
    double *post_h=NUMERIC_POINTER(R_post_h);
    double *post_yhat=NUMERIC_POINTER(R_post_yhat);
    //initialize to 0
    for(i=0;i<ng;i++){post_g[i]=0;post_b[i]=0;}
    for(i=0;i<nh;i++){post_h[i]=0;}
    for(i=0;i<n;i++){post_yhat[i]=0;}
    
   
 /*  //parameter to the power 2: needed for calculating SD.
    double post_mu2=0,post_var_e2=0,post_var_g2=0,post_var_b2=0,post_var_h2=0;
    double *post_b2= (double *) R_alloc(ng,sizeof(double));
    double *post_h2= (double *) R_alloc(nh,sizeof(double));
    double *post_yhat2=(double *) R_alloc(n,sizeof(double));
    for(i=0;i<ng;i++){
    post_b2[i]=0;
    }
    for(i=0;j<nh,i++){
    post_h2[i]=0;
    }
    for(i=0;i<n;i++){
    post_yhat2[i]=0;
    }
*/ 
    
    
   
    
    double *e=(double *) R_alloc(n,sizeof(double));
    double *X=(double *) R_alloc(n,sizeof(double));
    // including covariance matrix for g and b
    // SD.g was only done when L is NA, othersie, SD.delta_g is saved.
    double *XL, *XLh,*delta_g,*delta_b,*post_delta_g,*Xkb,*Xkg;
    
    if(!ISNAN(L[0])){
        //XL is the incidence matrix for delta_g
        XL=(double *) R_alloc(n*ng,sizeof(double));
        //XLh is the incidence matrix for delta_b
        XLh=(double *) R_alloc(n*ng,sizeof(double));
        delta_g=(double *)R_alloc(ng,sizeof(double));
        delta_b=(double *)R_alloc(ng,sizeof(double));
        post_delta_g=(double *)R_alloc(ng,sizeof(double));
      /*for SD.g
       // SEXP R_post_delta_g;
        //PROTECT(R_post_delta_g=allocVector(REALSXP,ng)); nProtect+=1;
        //post_delta_g=NUMERIC_POINTER(R_post_delta_g);
        //post_delta_g2=(double *)R_alloc(ng,sizeof(double));
        */
        for(j=0;j<ng;j++){
            post_delta_g[j]=0;
           // post_delta_g2[j]=0;
            for(i=0;i<n;i++){
            //XL is the incidence matrix for delta_g
            XL[n*j+i]=L[j*ng+C_IDL[i]];
            }
        }
        //Ldelta calculates b=L%*%delta_b or delta_b=Linv%*%b
        Ldelta(delta_g,Linv,g,ng);
        Ldelta(delta_b,Linv,b,ng);
    }
    //for SD.g
    //if(ISNAN(L[0])){
    //SEXP R_SD.g;
    //PROTECT(R_SD.g=allocVector(REALSXP,ng);nProtect+=1;
    //post_g2=(double *)R_alloc(ng,sizeof(double));
    //for(j=0;j<ng;j++){
    //post_g2[j]=0;
    //}
    //}
    
    
    //*initial values for e.//yStar is y except for the NA values;
    for(j=0;j<n;j++) e[j]=yStar[j]-mu[0]-g[C_IDL[j]]-(1+b[C_IDL[j]])*(h[C_IDE[j]]);
    /************************************************
     * //begin Gibbs sampler
     ************************************************/
    
    GetRNGstate();
    int sampleCount=0;
    for(i=0; i<nIter;i++){
        
        //sample environment effect h before b and g
        for(j=0;j<n;j++)X[j]=(1.0+b[C_IDL[j]]);
        sample_beta_ID(h,e,C_IDE,X,n,nh,var_e,var_h);
        
        double tXX;
        double tXy;

        if(ISNAN(L[0])){
            
            //sample b and g separately;
            //sample  b
            //	for(j=0;j<n;j++) X[j]=h[C_IDE[j]];
            //	sample_beta_ID(b,e,C_IDL,X,n,ng,var_e,var_b);
            //sample  g
            //	sample_beta_ID_x1(g,e,C_IDL,n,ng,var_e,var_g);
            
            //sample b and g together
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
            
        }else{
            //update XLh (incidence matrix for delta_b)
            for(j=0;j<n;j++) {
                for(k=0;k<ng;k++){
                    XLh[k*n+j]=XL[k*n+j]*h[C_IDE[j]];
                }
            }
            
            /*
             //sample b and g separately with L
             //sample delta_b
             sample_betaX(delta_b,e,XLh,n,ng,var_e,var_b);
             //sample delta_g
             sample_betaX(delta_g,e, XL, n, ng,var_e,var_g);
             */
            //sample b and g together
            
            for(k=0;k<ng;k++){
                tXy=0;
                tXX=0;
                Xkg=XL+k*n;
                Xkb=XLh+k*n;
                //sample b
                for(j=0;j<n;j++){
                    
                    e[j]=e[j]+delta_b[k]*Xkb[j];
                    tXX+=Xkb[j]*Xkb[j];
                    tXy+=Xkb[j]*e[j];
                    
                }
                delta_b[k]=sample_betaj(tXX,tXy,var_e,var_b);
                for(j=0;j<n;j++){
                    e[j]=e[j]-delta_b[k]*Xkb[j];
                }
                // sample g
                tXy=0;
                tXX=0;
                for(j=0;j<n;j++){
                    
                    e[j]=e[j]+delta_g[k]*Xkg[j];
                    tXX+=Xkg[j]*Xkg[j];
                    tXy+=Xkg[j]*e[j];
                    
                }
                delta_g[k]=sample_betaj(tXX,tXy,var_e,var_g);
                for(j=0;j<n;j++){
                    e[j]=e[j]-delta_g[k]*Xkg[j];
                }
            }
            
            //update b from delta_b
            Ldelta(b,L,delta_b,ng);
            //update g from delta_g //this is needed to get SD.g, but it is better not to do this for many lines, we can report delta_g
            //Ldelta(g,L,delta_g,ng);//
            
        }
        
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
        if(ISNAN(L[0])){
            SS=Sb;
            for(j=0;j<ng;j++) SS+=b[j]*b[j];
            DF=ng+dfb;
            var_b=SS/rchisq(DF);
            //var_g;
            SS=Sg;
            for(j=0;j<ng;j++)SS+=g[j]*g[j];
            DF=ng+dfg;
            var_g=SS/rchisq(DF);
        }else{
            for(j=0;j<ng;j++) SS+=delta_b[j]*delta_b[j];
            DF=ng+dfb;
            var_b=SS/rchisq(DF);
            //var_g;
            SS=Sg;
            for(j=0;j<ng;j++)SS+=delta_g[j]*delta_g[j];
            DF=ng+dfg;
            var_g=SS/rchisq(DF);
        }
        
        
        //sample intercept
        sample_mu(mu,e,var_e,n);
        
        //yhat & missing values
        for(j=0;j<n;j++){
            yhat[j]=yStar[j]-e[j];
        }
        if(nNa>0){
            for(j=0;j<nNa;j++){
                e[(whNA[j]-1)]=sqrtf(var_e)*norm_rand();
                yStar[(whNA[j]-1)]=yhat[(whNA[j]-1)]+e[(whNA[j]-1)];
            }
        }
        
        //running means and storing samples
        if(i>=(burnIn)){
            if((i+1)%thin==0){
            sampleCount+=1;
            post_mu += mu[0]/nSamples;
            post_var_e += var_e/nSamples;
            post_var_g += var_g/nSamples;
            post_var_b += var_b/nSamples;
            post_var_h += var_h/nSamples;
           /*//SD 
            post_mu2 += pow(mu[0],2)/nSamples;
            post_var_b2 += pow(var_b,2)/nSamples;
            post_var_e2 += pow(var_e,2)/nSamples;
            post_var_g2 += pow(var_g,2)/nSamples;
            post_var_h2 += pow(var_h,2)/nSamples;
            */
            
           
            
            for(j=0;j<nh;j++) {
            post_h[j] += h[j]/nSamples;
            //post_h2[j] += pow(h[j],2)/nSamples;
            }
            for(j=0;j<ng;j++) {
            post_b[j] += b[j]/nSamples;
            //post_b2[j] += pow(b[j],2)/nSamples;
            }
            
            //post_g
            if(ISNAN(L[0])){
                for(j=0;j<ng;j++){
                    post_g[j] += g[j]/nSamples;
                    //post_g2[j] += pow(g[j],2)/nSamples;
                }
            }else{
                for(j=0;j<ng;j++){
                    post_delta_g[j] += delta_g[j]/nSamples;
                    //post_delta_g2[j] += pow(delta_g[j],2)/nSamples;
                }
            }
           //post_yhat
           for(j=0;j<n;j++){
           post_yhat[j]+=yhat[j]/nSamples;
          // post_yhat2[j]+=pow(yhat[j],2)/nSamples;
           }
         //post_logLik
    	logLik=0;
    	if (nNa > 0) {
        for(j=0;j<nNotNa;j++){
        logLik += dnorm4(e[(whNotNA[j]-1)], 0, sqrtf(var_e), 1);
         }
        }else{
        for(j=0;j<n;j++){
        logLik += dnorm4(e[j], 0, sqrtf(var_e), 1);
        }
        }
        post_logLik += logLik/nSamples;


         //store samples in file
         
            fprintf(fsaveFile,"%f,%f,%f,%f,%f",mu[0],var_g,var_b,var_h,var_e);
            
            if(ISNAN(L[0])){
            	for(j=0;j<nVAR_Store;j++){
            	fprintf(fsaveFile,",%f",g[(VARstore[j]-1)]);
            	}
            }else{
            	for(j=0;j<nVAR_Store;j++){
            	fprintf(fsaveFile,",%f",delta_g[(VARstore[j]-1)]);
            	}
            }
            for(j=0;j<nVAR_Store;j++){
            	fprintf(fsaveFile,",%f",b[(VARstore[j]-1)]);
            }
            for(j=0;j<nENV_Store;j++){
            	fprintf(fsaveFile,",%f",h[(ENVstore[j]-1)]);
            }
            fprintf(fsaveFile,"\n");
    	  }
        
        }//end of running means and storing samples.
      
     
            // print out iterations
        if((i+1)%1000==0){Rprintf("iter:%d\n",i+1);}
        
    }//end iteration
    
    //printout the number saved samples and number of samples expected to save
   // Rprintf("nSamples:%d, SampleCount:%d\n",nSamples,sampleCount);
    //get post_g from post_delta_g
    if(!ISNAN(L[0])){Ldelta(post_g,L,post_delta_g,ng);}

    fclose(fsaveFile);

    PutRNGstate();//This write random numbers into R environments.

//return value to R

    REAL(R_post_mu)[0]=post_mu;
    REAL(R_post_var_g)[0]=post_var_g;
    REAL(R_post_var_b)[0]=post_var_b;
    REAL(R_post_var_h)[0]=post_var_h;
    REAL(R_post_var_e)[0]=post_var_e;
    //logLik
    REAL(R_post_logLik)[0]=post_logLik;
    double logLikAtPostMean;
    SEXP R_logLikAtPostMean;
    PROTECT(R_logLikAtPostMean=allocVector(REALSXP,1));nProtect+=1;
    double tmpE;
   logLikAtPostMean=0;
    if(nNa>0){
    for(j=0;j<nNotNa;j++){	
    tmpE=y[(whNotNA[j]-1)]-post_yhat[(whNotNA[j]-1)];
    logLikAtPostMean += dnorm4(tmpE, 0, sqrtf(post_var_e), 1);
    }
    }else{
    
    for(j=0;j<n;j++){	
    tmpE=y[j]-post_yhat[j];
    logLikAtPostMean += dnorm4(tmpE, 0, sqrtf(post_var_e), 1);
    }
    }
    
    REAL(R_logLikAtPostMean)[0]=logLikAtPostMean;
/*    //return standard deviations;
    SEXP R_SD.mu;
    PROTECT(R_SD.mu=allocVector(REALSXP,1);nProtect+=1;
    SEXP R_SD.var_e;
    PROTECT(R_SD.var_e=allocVector(REALSXP,1);nProtect+=1;
    SEXP R_SD.var_g;
    PROTECT(R_SD.var_g=allocVector(REALSXP,1);nProtect+=1;
    SEXP R_SD.var_b;
    PROTECT(R_SD.var_b=allocVector(REALSXP,1);nProtect+=1;
    SEXP R_SD.var_h;
    PROTECT(R_SD.var_h=allocVector(REALSXP,1);nProtect+=1;
    SEXP R_SD.b;
    PROTECT(R_SD.b=allocVector(REALSXP,ng);nProtect+=1;
    SEXP R_SD.h;
    PROTECT(R_SD.h=allocVector(REALSXP,nh);nProtect+=1;
    SEXP R_SD.yhat;
    PROTECT(R_SD.yhat=allocVector(REALSXP,n);nProtect+=1;
    REAL(R_SD.mu)[0]=sqrt(post_mu2-post_mu*post_mu);
    REAL(R_SD.var_g)[0]=sqrtf(post_var_g2-post_var_g*post_var_g);
    REAL(R_SD.var_b)[0]=sqrtf(post_var_b2-post_var_b*post_var_b);
    REAL(R_SD.var_h)[0]=sqrtf(post_var_h2-post_var_h*post_var_h);
    REAL(R_SD.var_e)[0]=sqrtf(post_var_e2-post_var_e*post_var_e);
    if(ISNAN(L[0])){
    for(j=0;j<ng;j++){
    REAL(R_SD.g)[j]=sqrtf(post_g2[j]-pow(post_g[j],2));
    }
    }else{
    for(j=0;j<ng;j++){
    REAL(R_SD.delta_g)[j]=sqrtf(post_delta_g2[j]-pow(post_delta_g[j],2));
    }
    }
    for(j=0;j<ng;j++){    
    REAL(R_SD.b)[j]=sqrtf(post_b2[j]-pow(post_b[j],2));
    }
    for(j=0;j<nh;j++){
    REAL(R_SD.h)[j]=sqrtf(post_h2[j]-pow(post_h[j],2));
    }
    for(j=0;j<n;j++){
    REAL(R_SD.yhat)[j]=sqrtf(post_yhat2[j]-pow(post_yhat[j],2));
    }
    // end of return SD
    */
 

    SEXP list;
    
    
    PROTECT(list = allocVector(VECSXP, 11));nProtect+=1;
   
    SET_VECTOR_ELT(list, 0, R_post_mu);
    SET_VECTOR_ELT(list, 1, R_post_var_g);
    SET_VECTOR_ELT(list, 2, R_post_var_b);
    SET_VECTOR_ELT(list, 3, R_post_var_h);
    SET_VECTOR_ELT(list, 4, R_post_var_e);
    SET_VECTOR_ELT(list, 5, R_post_g);
    SET_VECTOR_ELT(list, 6, R_post_b);
    SET_VECTOR_ELT(list,7, R_post_h);
    SET_VECTOR_ELT(list,8, R_post_yhat);
    SET_VECTOR_ELT(list,9,R_post_logLik);
    SET_VECTOR_ELT(list,10,R_logLikAtPostMean);
   /*/return SD. 
    SET_VECTOR_ELT(list,11,R_SD.mu);
    SET_VECTOR_ELT(list,12,R_SD.var_g);
    SET_VECTOR_ELT(list,13,R_SD.var_b);
    SET_VECTOR_ELT(list,14,R_SD.var_h);
    SET_VECTOR_ELT(list,15,R_SD.var_e);
    SET_VECTOR_ELT(list,16,R_SD.yhat);
    SET_VECTOR_ELT(list,17,R_SD.b);
    SET_VECTOR_ELT(list,18,R_SD.h);

    if(ISNAN(L[0])){
    SET_VECTOR_ELT(list,19,R_SD.g);
    }else{
    SET_VECTOR_ELT(list,20,R_post_delta_g);
    SET_VECTOR_ELT(list,21,R_SD.delta_g);
    }
   */ 

    UNPROTECT(nProtect);
   

    return(list);

}




