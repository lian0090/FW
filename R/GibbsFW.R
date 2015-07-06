GibbsFW=function(y,VAR,ENV,VARlevels=NULL,ENVlevels=NULL,savedir=".",nIter=5000,burnIn=3000,thin=1,df=5,dfg=5,dfh=5,dfb=5,S=NULL,Sg=NULL,Sb=NULL,Sh=NULL,A=NULL,inits=NULL,nchain=1,seed=NULL){
#check thin and df: they are functions in coda
  if(!is.numeric(thin)){
  	stop("thin must be a numeric")
  	}
  if(!is.numeric(df)){
  	stop("df must be a numeric")
  	}	
   current.dir=getwd()  
  if(!file.exists(savedir)){dir.create(savedir)}	
  setwd(savedir)
  ############################################# 
  # initialize
  ########################################################################################## 
  IDEL=getIDEL(VAR,ENV,VARlevels,ENVlevels)
  IDE=IDEL$IDE
  IDL=IDEL$IDL
  VARlevels=IDEL$VARlevels
  ENVlevels=IDEL$ENVlevels
  ng=length(VARlevels)
  nh=length(ENVlevels)
  inits=initialize(y,ng=ng,nh=nh,model="Gibbs",inits=inits,seed=seed,nchain=nchain)

  #hyper parameters:
  var_y=var(y)
  if(is.null(S))  S<-0.5*var_y*(df+2)  #S is the scale times df
  if(is.null(Sg)) Sg<-0.25*var_y*(dfg+2)
  if(is.null(Sb)) Sb<-0.5*sqrt(var_y)*(dfb+2)   
  if(is.null(Sh)) Sh<-0.5*sqrt(var_y)*(dfh+2)  
  if(!is.null(A)){
  	A=A[VARlevels,VARlevels]
  	L<-t(chol(A));
  	Linv=solve(L);
  	}else {
  		L<-NA;
  		Linv=NA;
  	}
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
      outi  =.Call("C_GibbsFW", y, IDL, IDE, g, b, h, nIter, burnIn, thin, sampFile,S, Sg, Sb, Sh, df, dfg, dfb, dfh, var_e, var_g, var_b, var_h, mu, as.vector(L), as.vector(Linv))
      names(outi)=c("mu","var_g","var_b","var_h","var_e","g","b","h");
      for(namej in names(outi)){
        assign(namej,outi[[namej]])
      }
      names(g)=VARlevels
      names(b)=VARlevels
      names(h)=ENVlevels
      postMean[[i]]=setFW(g=g,b=b,h=h,y=y,VAR=VAR,ENV=ENV,IDL=IDL,IDE=IDE,VARlevels=VARlevels,ENVlevels=ENVlevels,mu=mu,var_g=var_g,var_h=var_h,var_b=var_b,var_e=var_e)
    }
    names(postMean)=paste("Init",c(1:nchain),sep="")
    for(i in 1:nchain){
      sampFile=paste("sampsChain",i,".txt",sep="");
      samps[[i]] = mcmc(read.table(sampFile,sep=",", stringsAsFactors=F, header=T, check.names=F), start=burnIn+1, end=nIter, thin=thin)
      file.remove(sampFile)
    }
    	
    samps=mcmc.list(samps);	
    if(save_samps==TRUE){save(samps,file="Gibbs_samps.rda")}
    #mpsrf=gelman.diag(samps)$mpsrf
    #return(list(postMean=postMean,mpsrf=mpsrf))
    #setback the original working directory
    setwd(current.dir)
   
    return(list(postMean=postMean))
  }
  
  ############################################# 
  #end of runSampler function.
  #############################################
  #select seeds that produced higher correlation between ENVmean and hhat

  postMean=runSampler(inits=inits,nchain=nchain,nIter=nIter,burnIn=burnIn,save_samps=T,seed=seed)$postMean 
  class(postMean)=c("postMean","list")
  #cat(postMean[[1]]$corENVmean,"\n")
  save(postMean,file=file.path(savedir,"postMean_Gibbs.rda"))	
  setwd(current.dir)
  return(postMean=postMean)
}


.onUnload<-function(libpath){library.dynam.unload("GibbsFW",libpath)}

