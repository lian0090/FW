
#this is the Gibbs sampler written in C. R is just a wrapper.
GibbsFW=function(y,VAR,ENV,savedir=".",nIter=1000,burnIn=500,thin=1,df=5,dfg=5,dfh=5,dfb=5,S=NULL,Sg=NULL,Sb=NULL,Sh=NULL,A=NULL,inits=NULL,nchain=1,seed=NULL){
  
  current.dir=getwd()	
  if(!file.exists(savedir)){dir.create(savedir)}	
  setwd(savedir)
  #hyper parameters:
  var_y=var(y)
  if(is.null(S)) S<-0.5*var_y*(df+2)  #S is the scale times df
  if(is.null(Sg))Sg<-0.25*var_y*(dfg+2)
  if(is.null(Sb))Sb<-0.5*sqrt(var_y)*(dfb+2)  #should be sqrt (var(y)) 
  if(is.null(Sh))Sh<-0.5*sqrt(var_y)*(dfh+2)  #should be sqrt (var(y))
  if(!is.null(A)){L<-t(chol(A)); Linv=solve(L)}else {L<-NA;Linv=NA}
  ############################################# 
  # initialize
  ########################################################################################## 
  VAR=as.character(VAR)
  ENV=as.character(ENV)
  VARlevels=unique(VAR)
  ENVlevels=unique(ENV)
  fVAR=factor(VAR,levels=VARlevels,ordered=T)
  fENV=factor(ENV,levels=ENVlevels,ordered=T)
  IDL=as.numeric(fVAR)
  IDE=as.numeric(fENV)
  ng=length(VARlevels)
  nh=length(ENVlevels)
  inits=initialize(y,ng=ng,nh=nh,model="Gibbs",inits=inits,jags.seed=seed,nchain=nchain)
  
  
  ############################################# 
  # runSampler: A private function to run multiple chains
  #############################################
  runSampler=function(inits,nchain,nIter,burnIn,save_samps=F,seed=NULL){
    samps=list()	
    postMean=list()	
    for(i in 1:nchain){
      sampFile=paste("sampsChain",i,".txt",sep="");
      if(file.exists(sampFile)) file.remove(sampFile) else file.create(sampFile)
      mu=inits[[i]]$mu
      g=inits[[i]]$g
      b=inits[[i]]$b
      h=inits[[i]]$h
      var_e=inits[[i]]$var_e
      var_g=inits[[i]]$var_g
      var_b=inits[[i]]$var_b
      var_h=inits[[i]]$var_h
      if(!is.null(seed)){set.seed(seed)}
      postMean[[i]]  =.Call("C_GibbsFW", y, IDL, IDE, g, b, h, nIter, burnIn, thin, sampFile,S, Sg, Sb, Sh, df, dfg, dfb, dfh, var_e, var_g, var_b, var_h, mu, as.vector(L), as.vector(Linv))
      names(postMean[[i]])=c("mu","var_g","var_b","var_h","var_e","g","b","h");
      names(postMean[[i]]$g)=VARlevels
      names(postMean[[i]]$b)=VARlevels
      names(postMean[[i]]$h)=ENVlevels
      postMean[[i]]$fitted.values=postMean[[i]]$g[IDL]+(1+postMean[[i]]$b[IDL])*postMean[[i]]$h[IDE]
      postMean[[i]]$VARlevels=VARlevels
      postMean[[i]]$ENVlevels=ENVlevels
      postMean[[i]]$IDL=IDL
      postMean[[i]]$IDE=IDE
      class(postMean[[i]])=c("FW","list")
    }
    names(postMean)=paste("Init",c(1:nchain),sep="")
    for(i in 1:nchain){
      sampFile=paste("sampsChain",i,".txt",sep="");
      samps[[i]] = mcmc(read.table(sampFile,sep=",", stringsAsFactors=F, header=T, check.names=F), start=burnIn+1, end=nIter, thin=thin)
      file.remove(sampFile)
    }	
    samps=mcmc.list(samps);	
    if(save_samps==TRUE){save(samps,file="samps.RData")}
    #mpsrf=gelman.diag(samps)$mpsrf
    #return(list(postMean=postMean,mpsrf=mpsrf))
    return(list(postMean=postMean))
  }
  
  ############################################# 
  #end of runSampler function.
  #############################################
  
  postMean=runSampler(inits=inits,nchain=nchain,nIter=nIter,burnIn=burnIn,save_samps=T,seed=seed)$postMean   
  class(postMean)=c("postMean","list")
  save(postMean,file=file.path(savedir,"postMean_Gibbs.rda"))	
  setwd(current.dir)
  return(postMean=postMean)
}

.onUnload<-function(libpath){library.dynam.unload("GibbsFW",libpath)}

