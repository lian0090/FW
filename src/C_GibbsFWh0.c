//this code is for model y=mu+h+bh+g+e, where h~N(0,H\sigma^2_h)
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "sample_beta.h"
#include <stdlib.h>
#include <string.h>


SEXP C_GibbsFWh0(SEXP R_y, SEXP R_IDL, SEXP R_IDE, SEXP R_g, SEXP R_b, SEXP R_h, SEXP R_nIter, SEXP R_burnIn, SEXP R_thin, SEXP R_saveFile, SEXP R_S, SEXP R_Sg, SEXP R_Sb, SEXP R_Sh, SEXP R_df, SEXP R_dfg, SEXP R_dfb, SEXP R_dfh,SEXP R_var_e, SEXP R_var_g, SEXP R_var_b, SEXP R_var_h,SEXP R_mu,SEXP R_LA, SEXP R_LH, SEXP R_whNA , SEXP R_whNotNA, SEXP R_VARstore, SEXP R_ENVstore, SEXP R_yhatstore, SEXP R_Linv, SEXP R_LHinv)
{
    int nProtect=0;
    //data and initial values
    PROTECT(R_y=AS_NUMERIC(R_y));  nProtect+=1;
    PROTECT(R_IDL=AS_INTEGER(R_IDL)); nProtect+=1;
    PROTECT(R_IDE=AS_INTEGER(R_IDE)); nProtect+=1;
    PROTECT(R_g=AS_NUMERIC(R_g)); nProtect+=1;
    PROTECT(R_b=AS_NUMERIC(R_b)); nProtect+=1;
    PROTECT(R_h=AS_NUMERIC(R_h)); nProtect+=1;
    PROTECT(R_LA=AS_NUMERIC(R_LA)); nProtect+=1;
    PROTECT(R_LH=AS_NUMERIC(R_LH)); nProtect+=1;
    PROTECT(R_whNA=AS_INTEGER(R_whNA)); nProtect+=1;
    PROTECT(R_VARstore=AS_INTEGER(R_VARstore)); nProtect+=1;
    PROTECT(R_ENVstore=AS_INTEGER(R_ENVstore)); nProtect+=1;
    PROTECT(R_yhatstore=AS_INTEGER(R_yhatstore));nProtect+=1;
    PROTECT(R_Linv=AS_NUMERIC(R_Linv)); nProtect+=1;
    PROTECT(R_LHinv=AS_NUMERIC(R_LHinv)); nProtect+=1;
    int i,j,k;
    int n= length(R_y);
    int ng=length(R_g);
    int nh=length(R_h);
    int nNa=length(R_whNA);
    int nVAR_Store=length(R_VARstore);
    int nENV_Store=length(R_ENVstore);
    int nyhat_Store=length(R_yhatstore);
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
    int *IDL=INTEGER_POINTER(R_IDL);
    int *IDE=INTEGER_POINTER(R_IDE);
    double *L=NUMERIC_POINTER(R_LA);
    double *LH=NUMERIC_POINTER(R_LH);
    int *VARstore=INTEGER_POINTER(R_VARstore);
    int *ENVstore=INTEGER_POINTER(R_ENVstore);
    int *yhatstore=INTEGER_POINTER(R_yhatstore);
    double *Linv=NUMERIC_POINTER(R_Linv);
    double *LHinv=NUMERIC_POINTER(R_LHinv);

    //int MethodChol=INTEGER_VALUE(R_MethodChol);
    int ndeltah;
    int ndeltag;
    if(!ISNAN(L[0])){
    SEXP RdimUA;
    PROTECT(RdimUA=getAttrib(R_LA,R_DimSymbol));nProtect+=1;
    if(isNull(RdimUA)){ndeltag=1;}else ndeltag=INTEGER(RdimUA)[1];
    }else ndeltag=ng;
    //printf("ndeltag%d\n",ndeltag);
    
    if(!ISNAN(LH[0])){
    SEXP RdimUH;
    PROTECT(RdimUH=getAttrib(R_LH,R_DimSymbol));nProtect+=1;
    if(isNull(RdimUH)){ndeltah=1;}else ndeltah=INTEGER(RdimUH)[1];
    }else ndeltah=nh;
    //printf("ndeltah:%d\n",ndeltah);
    
    
    
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
    
    
    fprintf(fsaveFile,"%s,%s,%s,%s,%s","mu","var_g","var_b","var_h","var_e");
    
    
    for(j=0;j<nENV_Store;j++)fprintf(fsaveFile,",h[%d]",ENVstore[j]);
    
    
    for(j=0;j<nVAR_Store;j++){
        fprintf(fsaveFile,",b[%d]",VARstore[j]);
        fprintf(fsaveFile,",g[%d]",VARstore[j]);
    }
    
    for(j=0;j<nyhat_Store;j++){
        fprintf(fsaveFile,",yhat[%d]",yhatstore[j]);
    }
    
    fprintf(fsaveFile,"\n");
    
    
    /************************************************
     * posteria and yhat storage
     ************************************************/
    SEXP R_post_mu,R_post_var_g,R_post_var_b,R_post_var_h,R_post_var_e,R_post_g,R_post_b,R_post_h, R_post_yhat;
    PROTECT(R_post_mu=allocVector(REALSXP,1)); nProtect+=1;
    PROTECT(R_post_var_g=allocVector(REALSXP,1)); nProtect+=1;
    PROTECT(R_post_var_b=allocVector(REALSXP,1));nProtect+=1;
    PROTECT(R_post_var_h=allocVector(REALSXP,1));nProtect+=1;
    PROTECT(R_post_var_e=allocVector(REALSXP,1));nProtect+=1;
    PROTECT(R_post_g=allocVector(REALSXP,ng));nProtect+=1;
    PROTECT(R_post_b=allocVector(REALSXP,ng));nProtect+=1;
    PROTECT(R_post_h=allocVector(REALSXP,nh));nProtect+=1;
    PROTECT(R_post_yhat=allocVector(REALSXP,n));nProtect+=1;
    
    double post_mu=0, post_var_e=0,post_var_g=0,post_var_b=0,post_var_h=0;
    double *post_g=NUMERIC_POINTER(R_post_g);
    double *post_b=NUMERIC_POINTER(R_post_b);
    double *post_h=NUMERIC_POINTER(R_post_h);
    double *post_yhat=NUMERIC_POINTER(R_post_yhat);
    //initialize to 0
    for(i=0;i<ng;i++){post_g[i]=0;post_b[i]=0;}
    for(i=0;i<nh;i++){post_h[i]=0;}
    for(i=0;i<n;i++){post_yhat[i]=0;}
    
    //parameter to the power 2: needed for calculating SD_
    double post_mu2=0,post_var_e2=0,post_var_g2=0,post_var_b2=0,post_var_h2=0;
    double *post_g2= (double *) R_alloc(ng,sizeof(double));
    double *post_b2= (double *) R_alloc(ng,sizeof(double));
    double *post_h2= (double *) R_alloc(nh,sizeof(double));
    double *post_yhat2=(double *) R_alloc(n,sizeof(double));
    for(i=0;i<ng;i++){
        post_b2[i]=0;
        post_g2[i]=0;
    }
    for(i=0;i<nh;i++){
        post_h2[i]=0;
    }
    for(i=0;i<n;i++){
        post_yhat2[i]=0;
    }
    
    
    
    
    double *e=(double *) R_alloc(n,sizeof(double));
    double *X=(double *) R_alloc(n,sizeof(double));
    // including covariance matrix for g , b ,h
    double *ZgL, *ZgLh,*delta_g, *delta_b, *Xkb,*Xkg;
    //double *post_delta_g;
    double *ZhL, *ZhLb, *delta_h, *Xkh;
    if(!ISNAN(LH[0])){
        //ZhL : multiplier for delta_h
        ZhL=(double *) R_alloc(n*ndeltah,sizeof(double));
        //ZhLb is the incidence matrix for delta_h
        ZhLb=(double *) R_alloc(n*ndeltah,sizeof(double));
        delta_h=(double *)R_alloc(ndeltah,sizeof(double));
        //inital values for delta_h transformed from h.
        Ldelta(delta_h,LHinv,h,nh);
       // Udelta(delta_h,LHinv,h,ndeltah,nh);
        
        for(j=0;j<ndeltah;j++){
            //initial values for delta_h  same as h
            //delta_h[j]= NUMERIC_POINTER(R_h)[j];
            for(i=0;i<n;i++){
                ZhL[i+n*j]=LH[C_IDE[i]+j*nh];//only k=C_IDE[i] Zh_ik!=0, ZhL=LH[k,j]
            }
        }
    }
    
    if(!ISNAN(L[0])){
        //ZgL is the incidence matrix for delta_g
        ZgL=(double *) R_alloc(n*ndeltag,sizeof(double));
        //ZgLh is the incidence matrix for delta_b
        ZgLh=(double *) R_alloc(n*ndeltag,sizeof(double));
        delta_g=(double *)R_alloc(ndeltag,sizeof(double));
        delta_b=(double *)R_alloc(ndeltag,sizeof(double));
        //post_delta_g
        /*
        post_delta_g=(double *)R_alloc(ng,sizeof(double));
        for(j=0;j<ng;j++){
            post_delta_g[j]=0;
        }
        */
        
        //initial values for delta_g and delta_b
        //Ldelta calculates b=L%*%delta_b or delta_b=Linv%*%b
        //Udelta do the same thing, without recognizing L as lower triangular.
        Ldelta(delta_g,Linv,g,ng);
        Ldelta(delta_b,Linv,b,ng);
        //Udelta(delta_g,Linv,g,ndeltag,ng);
        //Udelta(delta_b,Linv,b,ndeltag,ng);
       
        for(j=0;j<ndeltag;j++){
            /*//initial values for delta_g and delta_b set the same as g and b
             delta_g[j]= NUMERIC_POINTER(R_g)[j];
             delta_b[j]= NUMERIC_POINTER(R_b)[j];
             */
            for(i=0;i<n;i++){
                //ZgL is the incidence matrix for delta_g
                ZgL[i+n*j]=L[C_IDL[i]+j*ng];
            }
        }
        
    }
    
    //*initial values for e.//yStar is y except for the NA values;
    for(j=0;j<n;j++) e[j]=yStar[j]-mu[0]-g[C_IDL[j]]-(1+b[C_IDL[j]])*(h[C_IDE[j]]);
    /************************************************
     * //begin Gibbs sampler
     ************************************************/
    double tXX;
    double tXy;
    GetRNGstate();
    int sampleCount=0;
    for(i=0; i<nIter;i++){
        
        //sample h
        if(ISNAN(LH[0])){
            for(j=0;j<n;j++)X[j]=(1.0+b[C_IDL[j]]);
            for(k=0;k<nh;k++){
                tXy=0;
                tXX=0;
                for(j=0;j<n;j++){
                    if(C_IDE[j]==k){
                        e[j]+=X[j]*h[k];
                        tXX+=pow(X[j],2);
                        tXy+=X[j]*e[j];
                    }
                }
                h[k]=sample_betaj(tXX,tXy,var_e,var_h);
            }
            for(j=0;j<n;j++){
                e[j]=e[j]-h[C_IDE[j]]*X[j];
            }
        }else{
            //update ZhLb (incidence matrix for delta_h)
            for(j=0;j<n;j++) {
                for(k=0;k<ndeltah;k++){
                    ZhLb[j+k*n]=ZhL[j+k*n]*(b[C_IDL[j]]+1);
                }
            }
            for(k=0;k<ndeltah;k++){
                tXy=0;
                tXX=0;
                Xkh=ZhLb+k*n;
                //sample b
                for(j=0;j<n;j++){
                    
                    e[j]=e[j]+delta_h[k]*Xkh[j];
                    tXX+=Xkh[j]*Xkh[j];
                    tXy+=Xkh[j]*e[j];
                    
                }
                delta_h[k]=sample_betaj(tXX,tXy,var_e,var_h);
                for(j=0;j<n;j++){
                    e[j]=e[j]-delta_h[k]*Xkh[j];
                }
            }
            //update h from delta_h
            Ldelta(h,LH,delta_h,nh);
           //Udelta(h,LH,delta_h,nh,ndeltah);
        }
        
        if(ISNAN(L[0])){
            
            //sample b and g
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
            //update ZgLh (incidence matrix for delta_b)
            for(j=0;j<n;j++) {
                for(k=0;k<ndeltag;k++){
                    ZgLh[k*n+j]=ZgL[k*n+j]*h[C_IDE[j]];
                }
            }
            for(k=0;k<ndeltag;k++){
                tXy=0;
                tXX=0;
                Xkg=ZgL+k*n;
                Xkb=ZgLh+k*n;
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
            //update g from delta_g //this is needed to get SD_g and plot the samples for g (instead of plot delta_g), but maybe it better not to do this for many lines for computation reaons, we might want to report delta_g in the case of too many line.
             Ldelta(b,L,delta_b,ng);
             Ldelta(g,L,delta_g,ng);
             //Udelta(g,L,delta_g,ng,ndeltag);
             //Udelta(b,L,delta_b,ng,ndeltag);
        }
        
        //var_e;
        SS=S;
        for(j=0;j<n;j++) SS+=e[j]*e[j];
        DF=n+df;
        var_e=SS/rchisq(DF);
        //var_h
        SS=Sh;
        if(ISNAN(LH[0])){
            for(j=0;j<nh;j++) SS+=h[j]*h[j];
            DF=nh+dfh;
            var_h=SS/rchisq(DF);
        }else{
            for(j=0;j<ndeltah;j++) SS+=pow(delta_h[j],2);
            DF=ndeltah+dfh;
            var_h=SS/rchisq(DF);
        }
        //var_b and var_g;
        
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
            SS=Sb;
            for(j=0;j<ndeltag;j++) SS+=delta_b[j]*delta_b[j];
            DF=ndeltag+dfb;
            var_b=SS/rchisq(DF);
            //var_g;
            SS=Sg;
            for(j=0;j<ng;j++)SS+=delta_g[j]*delta_g[j];
            DF=ndeltag+dfg;
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
                
                //SD
                post_mu2 += pow(mu[0],2)/nSamples;
                post_var_b2 += pow(var_b,2)/nSamples;
                post_var_e2 += pow(var_e,2)/nSamples;
                post_var_g2 += pow(var_g,2)/nSamples;
                post_var_h2 += pow(var_h,2)/nSamples;
                
                
                
                //post_h
                for(j=0;j<nh;j++) {
                    post_h[j] += h[j]/nSamples;
                    post_h2[j] += pow(h[j],2)/nSamples;
                }
                
                //post_b and post_g
                for(j=0;j<ng;j++){
                    post_g[j] += g[j]/nSamples;
                    post_g2[j] += pow(g[j],2)/nSamples;
                    post_b[j] += b[j]/nSamples;
                    post_b2[j] += pow(b[j],2)/nSamples;
                    
                }
                //post_yhat
                for(j=0;j<n;j++){
                    post_yhat[j]+=yhat[j]/nSamples;
                    post_yhat2[j]+=pow(yhat[j],2)/nSamples;
                }
                
                
                //store samples in file
                
                fprintf(fsaveFile,"%f,%f,%f,%f,%f",mu[0],var_g,var_b,var_h,var_e);
                for(j=0;j<nENV_Store;j++)fprintf(fsaveFile,",%f",h[(ENVstore[j]-1)]);
                for(j=0;j<nVAR_Store;j++){
                    fprintf(fsaveFile,",%f",b[(VARstore[j]-1)]);
                    fprintf(fsaveFile,",%f",g[(VARstore[j]-1)]);
                }
                for(j=0;j<nyhat_Store;j++){
                    fprintf(fsaveFile,",%f",yhat[(yhatstore[j]-1)]);
                }
                fprintf(fsaveFile,"\n");
            }
            
        }//end of running means and storing samples.
        
        
        // print out iterations
        if((i+1)%1000==0){Rprintf("iter:%d\n",i+1);}
        
    }//end iteration
    
    
    //get post_g from post_delta_g
   // if(!ISNAN(L[0])){
     //   Ldelta(post_g,L,post_delta_g,ng);
    //}
    
    fclose(fsaveFile);
    
    PutRNGstate();//This write random numbers into R environments.
    
    //return value to R
    
    REAL(R_post_mu)[0]=post_mu;
    REAL(R_post_var_g)[0]=post_var_g;
    REAL(R_post_var_b)[0]=post_var_b;
    REAL(R_post_var_h)[0]=post_var_h;
    REAL(R_post_var_e)[0]=post_var_e;
    
    //return standard deviations;
    SEXP R_SD_mu;
    PROTECT(R_SD_mu=allocVector(REALSXP,1));nProtect+=1;
    SEXP R_SD_var_e;
    PROTECT(R_SD_var_e=allocVector(REALSXP,1));nProtect+=1;
    SEXP R_SD_var_g;
    PROTECT(R_SD_var_g=allocVector(REALSXP,1));nProtect+=1;
    SEXP R_SD_var_b;
    PROTECT(R_SD_var_b=allocVector(REALSXP,1));nProtect+=1;
    SEXP R_SD_var_h;
    PROTECT(R_SD_var_h=allocVector(REALSXP,1));nProtect+=1;
    SEXP R_SD_g;
    PROTECT(R_SD_g=allocVector(REALSXP,ng));nProtect+=1;
    SEXP R_SD_b;
    PROTECT(R_SD_b=allocVector(REALSXP,ng));nProtect+=1;
    SEXP R_SD_h;
    PROTECT(R_SD_h=allocVector(REALSXP,nh));nProtect+=1;
    SEXP R_SD_yhat;
    PROTECT(R_SD_yhat=allocVector(REALSXP,n));nProtect+=1;
    REAL(R_SD_mu)[0]=sqrt(post_mu2-post_mu*post_mu);
    REAL(R_SD_var_g)[0]=sqrtf(post_var_g2-post_var_g*post_var_g);
    REAL(R_SD_var_b)[0]=sqrtf(post_var_b2-post_var_b*post_var_b);
    REAL(R_SD_var_h)[0]=sqrtf(post_var_h2-post_var_h*post_var_h);
    REAL(R_SD_var_e)[0]=sqrtf(post_var_e2-post_var_e*post_var_e);
    for(j=0;j<ng;j++){
        REAL(R_SD_g)[j]=sqrtf(post_g2[j]-pow(post_g[j],2));
        REAL(R_SD_b)[j]=sqrtf(post_b2[j]-pow(post_b[j],2));
    }
    
    for(j=0;j<nh;j++){
        REAL(R_SD_h)[j]=sqrtf(post_h2[j]-pow(post_h[j],2));
    }
    for(j=0;j<n;j++){
        REAL(R_SD_yhat)[j]=sqrtf(post_yhat2[j]-pow(post_yhat[j],2));
    }
    // end of return SD
    
    
    
    SEXP list;
    int countOut=0;
    int totalOut=18;
    
    PROTECT(list = allocVector(VECSXP, totalOut));nProtect+=1;
    
    SET_VECTOR_ELT(list, countOut, R_post_mu); countOut+=1;
    SET_VECTOR_ELT(list, countOut, R_post_var_g);countOut+=1;
    SET_VECTOR_ELT(list, countOut, R_post_var_b);countOut+=1;
    SET_VECTOR_ELT(list, countOut, R_post_var_h);countOut+=1;
    SET_VECTOR_ELT(list, countOut, R_post_var_e);countOut+=1;
    SET_VECTOR_ELT(list, countOut, R_post_g);countOut+=1;
    SET_VECTOR_ELT(list, countOut, R_post_b);countOut+=1;
    SET_VECTOR_ELT(list,countOut, R_post_h);countOut+=1;
    SET_VECTOR_ELT(list,countOut, R_post_yhat);countOut+=1;
    //return SD_
    SET_VECTOR_ELT(list,countOut,R_SD_mu);countOut+=1;
    SET_VECTOR_ELT(list,countOut,R_SD_var_g);countOut+=1;
    SET_VECTOR_ELT(list,countOut,R_SD_var_b);countOut+=1;
    SET_VECTOR_ELT(list,countOut,R_SD_var_h);countOut+=1;
    SET_VECTOR_ELT(list,countOut,R_SD_var_e);countOut+=1;
    SET_VECTOR_ELT(list,countOut,R_SD_g);countOut+=1;
    SET_VECTOR_ELT(list,countOut,R_SD_b);countOut+=1;
    SET_VECTOR_ELT(list,countOut,R_SD_h);countOut+=1;
    SET_VECTOR_ELT(list,countOut,R_SD_yhat);countOut+=1;
    
    
    UNPROTECT(nProtect);
    
    
    return(list);
    
}
