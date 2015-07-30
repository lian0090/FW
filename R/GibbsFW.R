GibbsFW=function(y,VAR,ENV,VARlevels=NULL,ENVlevels=NULL,saveAt=NULL,nIter=5000,burnIn=3000,thin=5,dfe=5,dfg=5,dfh=5,dfb=5,priorVARe=NULL,priorVARg=NULL,priorVARb=NULL,priorVARh=NULL,A=NULL,nchain=1,seed=NULL,inits=NULL,saveVAR=c(1:2),saveENV=c(1:2)){
#check Input type
	if(any(!is.numeric(c(nIter,burnIn,thin,dfe,dfg,dfh,dfb)))){
  		stop("thin and df must be a numeric")
  	}
	if(any(!is.numeric(c(saveVAR,saveENV)))){
  		stop("saveVAR and saveENV must be a numeric")
  	}	
	if(is.null(saveAt)){
		saveAt=paste(getwd(),"/",sep="")
	} 
	savedir=gsub("/[^/]*$","",saveAt)
   
	if(!file.exists(savedir)){
		dir.create(savedir)
	}	
  
	if(!is.null(seed)){
		if(length(seed)!=nchain)stop("number of seed must be equal to the number of chains")
	}	
  
	if(!is.integer(saveVAR)){
		stop("saveVAR must be interger vector")
  	}
	if(!is.integer(saveENV)){
		stop("saveENV must be interger vector")
  	}	
	IDEL=getIDEL(VAR,ENV,VARlevels,ENVlevels)
	IDE=IDEL$IDE
	IDL=IDEL$IDL
	VARlevels=IDEL$VARlevels
	ENVlevels=IDEL$ENVlevels
  	
  ############################################# 
  # initialize
  ########################################################################################## 
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
  	 	
  	 	if(any(!(unique(VAR) %in% colnames(A)))) stop("Covariance structure not available for some varieties")
  	 
  		if(is.null(VARlevels)){
  	 		VARlevels=colnames(A)  	 
  	 	}else{
  			A=A[VARlevels,VARlevels] ##changed the order of A to correspond to user defined or default VARlevels.
  		}	
  	}	
 
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
	if(!is.null(A)){
		L<-t(chol(A));
		#Linv=solve(L);
		Linv=forwardsolve(L,x=diag(1,nrow(L)),upper.tri=F)
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
  	
  		#post_logLik
  		#postMeanlogLikT=logLikAtPostMeanT=rep(NA,nchain)
  	
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
			if(!is.null(seed)){set.seed(seed[i])}
			outi  =.Call("C_GibbsFW", y, IDL, IDE, g, b, h, nIter, burnIn, thin, sampFile,S, Sg, Sb, Sh, dfe, dfg, dfb, dfh, var_e, var_g, var_b, var_h, mu, as.vector(L), as.vector(Linv),whNA,whNotNA,saveVAR,saveENV)
			names(outi)=c("mu","var_g","var_b","var_h","var_e","g","b","h","post_yhat");
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
			# postMeanlogLikT[i]=outi$postlogLik
			# logLikAtPostMeanT[i]=outi$logLikAtPostMean
	
    	}
    
		rownames(gT)=VARlevels
		rownames(bT)=VARlevels
		rownames(hT)=ENVlevels
		colnames(post_yhatT) = colnames(gT) = colnames(bT) = colnames(hT) = names(muT) = names(var_eT) = names(var_gT) = names(var_bT) = names(var_hT) = paste("Init",c(1:nchain),sep="")

 		yhatT=muT+gT[VAR,,drop=F]+(1+bT[VAR,,drop=F])*hT[ENV,,drop=F] 
      
  
		#but DIC is not a good indicator of prediction accuracy. it seesm that standardizing by saturated likelihood (when the variance unkown) helps.But we are not going to test this for now.  
		#names(postMeanlogLikT)=names(logLikAtPostMeanT)= paste("Init",c(1:nchain),sep="")   
		# pD=-2*(postMeanlogLikT-logLikAtPostMeanT)
		# DIC=pD-2*postMeanlogLikT
		# fit=list(postMeanlogLik=postMeanlogLikT,logLikAtPostMean=logLikAtPostMeanT,pD=pD,DIC=DIC)
		# postMean=list(y=y,whichNa=whNA,VAR=VAR,ENV=ENV,VARlevels=VARlevels,ENVlevels=ENVlevels, mu=muT,g=gT,b=bT,h=hT,yhat=yhatT,var_e=var_eT,var_g=var_gT,var_b=var_bT,var_h=var_hT,post_yhat=post_yhatT,fit=fit)
    
		 postMean=list(y=y,whichNa=whNA,VAR=VAR,ENV=ENV,VARlevels=VARlevels,ENVlevels=ENVlevels, mu=muT,g=gT,b=bT,h=hT,yhat=yhatT,var_e=var_eT,var_g=var_gT,var_b=var_bT,var_h=var_hT,post_yhat=post_yhatT)
    
		class(postMean)=c("FW","list")
    
		for(i in 1:nchain){
			sampFile=paste("sampsChain",i,".txt",sep="");
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

