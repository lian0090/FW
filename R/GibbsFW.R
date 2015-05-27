#this is the Gibbs sampler written in C. R is just a wrapper.
GibbsFW=function(y,IDL,IDE,savedir=".",nIter=1000,burnIn=500,thin=1,df=5,dfg=5,dfh=5,dfb=5,S=NULL,Sg=NULL,Sb=NULL,Sh=NULL,A=NULL,inits=NULL,nchain=1,seed=NULL){
            		
	current.dir=getwd()	
	if(!file.exists(savedir)){dir.create(savedir)}	
	setwd(savedir)
	#hyper parameters:
	if(is.null(S)) S<-0.5*var(y)*(df+2)  #S is the scale times df
	if(is.null(Sg))Sg<-0.25*var(y)*(dfg+2)
	if(is.null(Sb))Sb<-0.5*var(y)*(dfb+2)  
	if(is.null(Sh))Sh<-0.5*var(y)*(dfh+2)
	if(!is.null(A)){L<-t(chol(A)); Linv=solve(L)}else {L<-NA;Linv=NA}
############################################# 
# initialize
########################################################################################## 
	if((!is.integer(IDE)) | (!is.integer(IDL))){stop("IDE and IDL must be integers")}
	ng=length(unique(IDL))
	nh=length(unique(IDE))
	inits=initialize(y,ng=ng,nh=nh,model="Gibbs",inits=inits,jags.seed=seed,nchain=nchain)
	
  
############################################# 
# Run Gibbs Sampler for each chain
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
		}
		names(postMean)=paste("Init",c(1:nchain),sep="")
		for(i in 1:nchain){
			sampFile=paste("sampsChain",i,".txt",sep="");
			samps[[i]] = mcmc(read.table(sampFile,sep=",", stringsAsFactors=F, header=T, check.names=F), start=burnIn+1, end=nIter, thin=thin)
			file.remove(sampFile)
		}	
		samps=mcmc.list(samps);	
		if(save_samps==TRUE){save(samps,file="samps.rda")}
		#mpsrf=gelman.diag(samps)$mpsrf
		#return(list(postMean=postMean,mpsrf=mpsrf))
		return(list(postMean=postMean))
}
##use multiple chains to choose the seed for best convergence.	    
   # mpsrf=vector();
   #newseed=vector()
#	if(selectSeed==T){
#		if(nchain>1){
# 			for(i in 1:10){
#    		seed=sample(1:1000,1) 
#    		newseed[i]=seed		
#    		mpsrf[i]=runSampler(inits=inits,nchain=nchain,nIter=1000,burnIn=500,save_samps=F,seed=seed)$mpsrf
#    		cat("testing seed=",seed,", mpsrf=",mpsrf[i],"\n")
#    	}
#    		seed=newseed[which.min(mpsrf)]	
#    	}else{warning("seed selection can only be done for more than 1 chain")}
#    	}
#if(ans$mpsrf>1.3){warning("Gelman and Rubin’s multivariate potential scale reduction factor is ", formatC(ans$mpsrf,digits=2,format="f"))}  

postMean=runSampler(inits=inits,nchain=nchain,nIter=nIter,burnIn=burnIn,save_samps=T,seed=seed)$postMean	 
save(postMean,file=file.path(savedir,"postMean_Gibbs.rda"))	

setwd(current.dir)
return(postMean=postMean)
}

.onUnload<-function(libpath){library.dynam.unload("GibbsFW",libpath)}

