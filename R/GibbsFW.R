GibbsFW=function(y,VAR,ENV,VARlevels=NULL,ENVlevels=NULL,savedir=".",nIter=5000,burnIn=3000,thin=1,df=5,dfg=5,dfh=5,dfb=5,S=NULL,Sg=NULL,Sb=NULL,Sh=NULL,A=NULL,inits=NULL,nchain=1,seed=NULL){
#check thin and df: they are functions in coda
  if(!is.numeric(thin)){
  	stop("thin must be a numeric")
  	}
  if(!is.numeric(df)){
  	stop("df must be a numeric")
  	}	
   current.dir=getwd()  
  if(!file.exists(savedir)){
  	dir.create(savedir)
  	}	
  setwd(savedir)
  ############################################# 
  # initialize
  ########################################################################################## 
  if(is.null(VARlevels)){
  	if(!is.null(A)){
  	 if(is.null(colnames(A))){
  	 	stop("Variety names for variance matrix must be provided as column names")
  	 }
  	 if(!is.null(VARlevels)){
  	 	stop("When covariance matrix A is provided, VARlevels are automatically the column names of A, user should not specify VARlevels in this case")
  	 }
  		VARlevels=colnames(A)
  	}	
  }
  IDEL=getIDEL(VAR,ENV,VARlevels,ENVlevels)
  IDE=IDEL$IDE
  IDL=IDEL$IDL
  VARlevels=IDEL$VARlevels
  ENVlevels=IDEL$ENVlevels
  ng=length(VARlevels)
  nh=length(ENVlevels)
  whNA=which(is.na(y))
  nNa=length(whNA)
  
  inits=initialize.Gibbs(y,ng=ng,nh=nh,inits=inits,seed=seed,nchain=nchain)

  #hyper parameters:
  var_y=var(y,na.rm=T)
  if(is.null(S))  S<-0.5*var_y*(df+2)  #S is the scale times df
  if(is.null(Sg)) Sg<-0.25*var_y*(dfg+2)
  if(is.null(Sb)) Sb<-0.5*sqrt(var_y)*(dfb+2)   
  if(is.null(Sh)) Sh<-0.5*sqrt(var_y)*(dfh+2)  
  if(!is.null(A)){
  	A=A[VARlevels,VARlevels] #when VARlevels was specifed by the user.
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
  	samps=vector("list",nchain)
  	gT=bT=matrix(NA,ng,nchain)
  	hT=matrix(NA,nh,nchain)
  	var_eT=var_gT=var_bT=var_hT=muT=rep(NA,nchain)
  	post_yhatT=matrix(NA,nrow=length(y),ncol=nchain)
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
      outi  =.Call("C_GibbsFW", y, IDL, IDE, g, b, h, nIter, burnIn, thin, sampFile,S, Sg, Sb, Sh, df, dfg, dfb, dfh, var_e, var_g, var_b, var_h, mu, as.vector(L), as.vector(Linv),whNA)
      names(outi)=c("mu","var_g","var_b","var_h","var_e","g","b","h","post_yhat");
      gT[,i]=outi$g
      bT[,i]=outi$b
      hT[,i]=outi$h
      var_eT[i]=outi$var_e
      var_gT[i]=outi$var_g
      var_bT[i]=outi$var_b
      var_hT[i]=outi$var_h
      muT[i]=outi$mu
      post_yhatT[,i]=outi$post_yhat
   
    }
    
      rownames(gT)=VARlevels
      rownames(bT)=VARlevels
      rownames(hT)=ENVlevels
      colnames(post_yhatT)=colnames(gT)=colnames(bT)=colnames(hT)=names(muT)=names(var_eT)=names(var_gT)=names(var_bT)=names(var_hT)=paste("Init",c(1:nchain),sep="")
  
  
     yhatT=muT+gT[VAR,]+(1+bT[VAR,])*hT[ENV,]
      
     postMean=list(y=y,whichNa=whNA,VAR=VAR,ENV=ENV,VARlevels=VARlevels,ENVlevels=ENVlevels, mu=muT,g=gT,b=bT,h=hT,yhat=yhatT,var_e=var_eT,var_g=var_gT,var_b=var_bT,var_h=var_hT,post_yhat=post_yhatT)
    class(postMean)=c("FW","list")
    
     for(i in 1:nchain){
      sampFile=paste("sampsChain",i,".txt",sep="");
      samps[[i]] = mcmc(read.table(file.path(savedir,sampFile),sep=",", stringsAsFactors=F, header=T, check.names=F), start=burnIn+1, end=nIter, thin=thin)
      file.remove(sampFile)
    }
    	
    samps=mcmc.list(samps);	
	if(save_samps==TRUE){
		save(samps,file="Gibbs_samps.rda")
		}
    
    #mpsrf=gelman.diag(samps)$mpsrf
    #return(list(postMean=postMean,mpsrf=mpsrf))
    #setback the original working directory
    setwd(current.dir)
   
    return(postMean)
  }
  
  ############################################# 
  #end of runSampler function.
  #############################################
  postMean=runSampler(inits=inits,nchain=nchain,nIter=nIter,burnIn=burnIn,save_samps=T,seed=seed)
  #save(postMean,file=file.path(savedir,"postMean_Gibbs.rda"))	
  setwd(current.dir)
  return(postMean)
}



.onUnload<-function(libpath){library.dynam.unload("GibbsFW",libpath)}

