#a wrapper for GibbsFW and lmFW
FW=function(y,VAR,ENV,method=c("OLS","Gibbs")[2], savedir=".",nIter=5000,burnIn=3000,thin=1,df=5,dfg=5,dfh=5,dfb=5,S=NULL,Sg=NULL,Sb=NULL,Sh=NULL,A=NULL,inits=NULL,nchain=1,seed=NULL){

whNA=which(is.na(y))
#if genotype or environment is completely missing for a GxE combination, the predicted value of  y is still NA.
  

if(method=="OLS"){
	 predictedValue=lmFW(y,VAR,ENV,savedir=savedir)
}

if(method=="Gibbs"){
	 	 predictedValue=GibbsFW(y=y,VAR=VAR,ENV=ENV,nIter=nIter,burnIn=burnIn,thin=thin,df=df,dfg=dfg,dfh=dfh,dfb=dfb,S=S,Sg=Sg,Sb=Sb,Sh=Sh, A=A,inits=inits,nchain=nchain,seed=seed,savedir=savedir);	
	
}

return(predictedValue)

}
