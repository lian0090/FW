GibbsFWh0=function(y,VAR,ENV,saveAt="",nIter=5000,burnIn=3000,thin=5,dfe=5,dfg=5,dfh=5,dfb=5,priorVARe=NULL,priorVARg=NULL,priorVARb=NULL,priorVARh=NULL,A=NULL,H=NULL,nchain=1,seed=NULL,inits=NULL,saveVAR=c(1:2),saveENV=c(1:2),saveyhat=NULL){
  #check Input type
  if(any(!is.numeric(c(nIter,burnIn,thin,dfe,dfg,dfh,dfb)))){
    stop("thin and df must be a numeric")
  }
  
  if (saveAt == "") {
    saveAt = paste(getwd(), "/", sep = "")
  }
  ##paste a randomfile name to saveAt in order to get the correct directory name.
  #for example, dirname("~/a//") will be ~. 
  savedir=dirname(paste(saveAt,"ArandomfileName.txt",sep=""))
  
  if(!file.exists(savedir)){
    dir.create(savedir)
  }	
  
  if(!is.null(seed)){
    if(length(seed)!=nchain)stop("number of seed must be equal to the number of chains")
  }	
  VAR=factor(VAR)
  ENV=factor(ENV)
  IDL=as.numeric(VAR)
  IDE=as.numeric(ENV)
  VARlevels=levels(VAR)
  ENVlevels=levels(ENV)
  if(!is.integer(saveVAR) & !(is.character(saveVAR))){
    stop("saveVAR must be interger or character vector\n")
    if(is.character(saveVAR))saveVAR=match(saveVAR,VARlevels,nomatch=0)
  
  }
  if(!is.integer(saveENV) & !(is.character(saveENV))){
    stop("saveENV must be interger or character vector\n")
    if(is.character(saveENV))saveENV=match(saveENV,ENVlevels,nomatch=0)
  }  
  
  if(is.null(saveyhat))saveyhat=which((IDL %in% saveVAR) & (IDE %in% saveENV))
  
  if(!is.null(H)){
    if(nrow(H)!=ncol(H)){
      stop("H must be square matrix")
    }
    if(is.null(colnames(H))| is.null(rownames(H))){
      stop("Environment names for variance matrix must be provided as column names and rownames\n")
    }
    
    if(!all(colnames(H)==rownames(H))){
      stop("colnames of H must be equal to rownames of H\n")
    }
    
    if(any(!(ENVlevels %in% colnames(H)))) stop("Covariance structure not available for some environments")
    
    H=H[ENVlevels,ENVlevels] ##changed the order of A to correspond to user defined or default VARlevels.
    
  }	
  
  
  
  if(!is.null(A)){
    if(nrow(A)!=ncol(A)){
      stop("A must be square matrix")
    }
    if(is.null(colnames(A))| is.null(rownames(A))){
      stop("Variety names for variance matrix must be provided as column names and rownames\n")
    }
    
    if(!all(colnames(A)==rownames(A))){
      stop("colnames of A must be equal to rownames of A\n")
    }
    
    if(any(!(VARlevels %in% colnames(A)))) stop("Covariance structure not available for some varieties")
    
    
    A=A[VARlevels,VARlevels] ##changed the order of A to correspond to user defined or default VARlevels.
    
  }	
  
  ############################################# 
  # initialize
  ########################################################################################## 
  
  
  ng=length(VARlevels)
  nh=length(ENVlevels)
  whNA=which(is.na(y))
  n=length(y)
  whNotNA=setdiff(1:n,whNA)
  nNa=length(whNA)
  if(max(saveVAR)>ng | min(saveVAR)<1){stop("saveVAR must be in the value of 1 to ng")}
  if(max(saveENV)>nh | min(saveENV)<1){stop("saveENV must be in the value of 1 to nh")}
  inits=initialize.Gibbs(y,ng=ng,nh=nh,inits=inits,nchain=nchain)
  
  #hyper parameters:
  var_y=var(y,na.rm=T)
  if(is.null(priorVARe)) {priorVARe=0.5*var_y}; S<-priorVARe*(dfe+2)  #S is the scale times df
  if(is.null(priorVARg)) {priorVARg=0.25*var_y};Sg<-priorVARg*(dfg+2)
  if(is.null(priorVARb)) {priorVARb=0.5*sqrt(var_y)}; Sb<-priorVARb*(dfb+2)   
  if(is.null(priorVARh)) {priorVARh=0.5*sqrt(var_y)}; Sh<-priorVARh*(dfh+2)  
  #if(!is.null(A)) {LA<-t(chol(A))} else {LA=NA}
  #if(!is.null(H)) {LH=t(chol(H))} else {LH=NA}
  ##LAinv and LHinv is only used to get initial values for delta_h, delta_b and delta_g
  if(!is.null(A)) {
    LA<-try(t(chol(A)),silent=T); 
    if(inherits(LA,"try-error")) stop("A is not full rank, try add small values to the diagonal of A to make it full rank.  ");
    LAinv=forwardsolve(LA,x=diag(1,nrow(LA)),upper.tri=F)
  } else {
    LA=NA;
    LAinv=NA
  }
  if(!is.null(H)) {
    LH<-try(t(chol(H)),silent=T);  
    if(inherits(LH,"try-error")) stop("H is not full rank, try add small values to the diagonal of H to make it full rank."); 
    LHinv=forwardsolve(LH,x=diag(1,nrow(LH)),upper.tri=F)
  } else {
    LH=NA;
    LHinv=NA
  }
  # this is used for eigen decomposition method.   
  # 	if(!is.null(A)) {
  #     eigenA=eigen(A,symmetric=T); whichEigen=which(eigenA$values>1e-7); 
  #     U1=eigenA$vectors[,whichEigen]
  #    d1=eigenA$values[whichEigen]
  #    LA=U1%*%diag(sqrt(d1))
  #    LAinv=diag(sqrt(1/d1))%*%t(U1)
  # 	}else {LA=NA;LAinv=NA}
  #  	 if(!is.null(H)) {
  #  	   eigenH=eigen(H,symmetric=T); whichEigen=which(eigenH$values>1e-7); 
  #  	   U1=eigenH$vectors[,whichEigen]
  #  	   d1=eigenH$values[whichEigen]
  # 	   LH=U1%*%diag(sqrt(d1))
  #  	   LHinv=diag(sqrt(1/d1))%*%t(U1)
  #  	 } else {LH=NA;LHinv=NA}
  ############################################# 
  # runSampler: A private function to run multiple chains
  #############################################
  runSampler=function(inits,nchain,nIter,burnIn,save_samps=F,seed=NULL){
    samps=vector("list",nchain)
    gT=bT=matrix(NA,ng,nchain)
    hT=matrix(NA,nh,nchain)
    var_eT=var_gT=var_bT=var_hT=muT=rep(NA,nchain)
    post_yhatT=matrix(NA,nrow=length(y),ncol=nchain)
    ##standard deviations
    SD.gT=SD.bT=matrix(NA,ng,nchain)
    SD.hT=matrix(NA,nh,nchain)
    SD.var_eT=SD.var_gT=SD.var_bT=SD.var_hT=SD.muT=rep(NA,nchain)
    SD.yhatT=matrix(NA,nrow=length(y),ncol=nchain)
    
    
    
    #post_logLik
    #postMeanlogLikT=logLikAtPostMeanT=rep(NA,nchain)
    
    for(i in 1:nchain){
      sampFile=paste(saveAt, "sampsChain",i,".txt",sep="");
      sampFile=normalizePath(sampFile,mustWork=F)
      file.create(sampFile)
      mu=inits[[i]]$mu
      g=inits[[i]]$g
      b=inits[[i]]$b
      h=inits[[i]]$h
      var_e=inits[[i]]$var_e
      var_g=inits[[i]]$var_g
      var_b=inits[[i]]$var_b
      var_h=inits[[i]]$var_h
      if(!is.null(seed)){set.seed(seed[i])}
      outi  =.Call("C_GibbsFWh0", y, IDL, IDE, g, b, h, nIter, burnIn, thin, sampFile,S, Sg, Sb, Sh, dfe, dfg, dfb, dfh, var_e, var_g, var_b, var_h, mu, LA, LH,whNA,whNotNA,saveVAR,saveENV,saveyhat,LAinv,LHinv)
      names(outi)=names(outi)=c("mu","var_g","var_b","var_h","var_e","g","b","h","post_yhat","SD.mu","SD.var_g","SD.var_b","SD.var_h","SD.var_e","SD.g","SD.b","SD.h","SD.yhat");
      #when there is postlogLik
      #names(outi)=c("mu","var_g","var_b","var_h","var_e","g","b","h","post_yhat","postlogLik","logLikAtPostMean");
      gT[,i]=outi$g
      bT[,i]=outi$b
      hT[,i]=outi$h
      var_eT[i]=outi$var_e
      var_gT[i]=outi$var_g
      var_bT[i]=outi$var_b
      var_hT[i]=outi$var_h
      muT[i]=outi$mu
      post_yhatT[,i]=outi$post_yhat
      ##standard deviations
      SD.gT[,i]=outi$SD.g
      SD.bT[,i]=outi$SD.b
      SD.hT[,i]=outi$SD.h
      SD.var_eT[i]=outi$SD.var_e
      SD.var_gT[i]=outi$SD.var_g
      SD.var_bT[i]=outi$SD.var_b
      SD.var_hT[i]=outi$SD.var_h
      SD.muT[i]=outi$SD.mu
      SD.yhatT[,i]=outi$SD.yhat
      
      
      # postMeanlogLikT[i]=outi$postlogLik
      # logLikAtPostMeanT[i]=outi$logLikAtPostMean
    }
    
    rownames(gT)=VARlevels
    rownames(bT)=VARlevels
    rownames(hT)=ENVlevels
    colnames(post_yhatT) = colnames(gT) = colnames(bT) = colnames(hT) = names(muT) = names(var_eT) = names(var_gT) = names(var_bT) = names(var_hT) = paste("Init",c(1:nchain),sep="")
    
    ##rownames and colnames for SD
    rownames(SD.gT)=VARlevels
    rownames(SD.bT)=VARlevels
    rownames(SD.hT)=ENVlevels
    colnames(SD.yhatT)=colnames(SD.gT)=colnames(SD.bT)=colnames(SD.hT)=names(SD.muT)=names(SD.var_eT)=names(SD.var_gT)=names(SD.var_bT)=names(SD.var_hT)=paste("Init",c(1:nchain),sep="")
    
    #yhatT=muT+gT[VAR,,drop=F]+(1+bT[VAR,,drop=F])*hT[ENV,,drop=F] 
    
    
    #but DIC is not a good indicator of prediction accuracy. it seesm that standardizing by saturated likelihood (when the variance unkown) helps.But we are not going to test this for now.  
    #names(postMeanlogLikT)=names(logLikAtPostMeanT)= paste("Init",c(1:nchain),sep="")   
    # pD=-2*(postMeanlogLikT-logLikAtPostMeanT)
    # DIC=pD-2*postMeanlogLikT
    # fit=list(postMeanlogLik=postMeanlogLikT,logLikAtPostMean=logLikAtPostMeanT,pD=pD,DIC=DIC)
    # postMean=list(y=y,whichNa=whNA,VAR=VAR,ENV=ENV,VARlevels=VARlevels,ENVlevels=ENVlevels, mu=muT,g=gT,b=bT,h=hT,yhat=yhatT,var_e=var_eT,var_g=var_gT,var_b=var_bT,var_h=var_hT,post_yhat=post_yhatT,fit=fit)
    
    postMean=list(y=y,whichNa=whNA,VAR=VAR,ENV=ENV, VARlevels=VARlevels,ENVlevels=ENVlevels, mu=muT,SD.mu=SD.muT,g=gT,SD.g=SD.gT,b=bT,SD.b=SD.bT,h=hT,SD.h=SD.hT,yhat=post_yhatT,SD.yhat=SD.yhatT,var_e=var_eT,SD.var_e=SD.var_eT,var_g=var_gT,SD.var_g=SD.var_gT,var_b=var_bT,SD.var_b=SD.var_bT,var_h=var_hT,SD.var_h=SD.var_hT)		
    class(postMean)=c("FW","list")
    
    for(i in 1:nchain){
      sampFile=paste(saveAt, "sampsChain",i,".txt",sep="");
      samps[[i]] = mcmc(read.table(sampFile,sep=",", stringsAsFactors=F, header=T, check.names=F), start=burnIn+thin, end=nIter, thin=thin)
      file.remove(sampFile)
    }
    
    samps=mcmc.list(samps);	
    if(save_samps==TRUE){
      save(samps,file=paste(saveAt,"samps.rda",sep=""))
    }
    
    #mpsrf=gelman.diag(samps)$mpsrf
    #return(list(postMean=postMean,mpsrf=mpsrf))
    
    return(postMean)
  }
  
  #end of runSampler function.
  
  postMean=runSampler(inits=inits,nchain=nchain,nIter=nIter,burnIn=burnIn,save_samps=T,seed=seed)
  #save(postMean,file=file.path(savedir,"postMean_Gibbs.rda"))	
  
  return(postMean)
}




.onUnload<-function(libpath){library.dynam.unload("FW",libpath)}

