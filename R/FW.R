#a wrapper for GibbsFW and lmFW
FW=function(y,VAR,ENV,method=c("OLS","Gibbs")[2], VARlevels=NULL,ENVlevels=NULL,savedir=".",nIter=5000,burnIn=3000,thin=1,df=5,dfg=5,dfh=5,dfb=5,S=NULL,Sg=NULL,Sb=NULL,Sh=NULL,A=NULL,inits=NULL,nchain=1,seed=NULL){

if(method=="OLS"){
	 predictedValue=lmFW(y,VAR,ENV,VARlevels,ENVlevels,savedir=savediri)
}

if(method=="Gibbs"){
	 	 predictedValue=GibbsFW(y=y,VAR=VAR,ENV=ENV,VARlevels=VARlevels,ENVlevels=ENVlevels,nIter=nIter,burnIn=burnIn,thin=thin,A=NULL,seed=seed,savedir=savediri);	
	
}

return(predictedValue)

}
