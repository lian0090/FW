#initialize
initialize.Gibbs=function(y,ng,nh,inits=NULL,nchain=1){
  whichNa=which(is.na(y))
  #currently not using, need to add input variable VAR, and ENV to use the following
  #ynNA=y[-whichNa]
  #VAR_nNA=VAR[-whichNa]
  #ENV_nNA=ENV[-whichNa]
  #ng_nNA=length(unique(VAR_nNA))
  #nh_nNA=length(unique(ENV_nNA))
  #seed is to set the random seed for Gibbs Sampler for jags. Not for the random seed of inital values.
  #the seed for setting up initial values is their chain index
  var_y=var(y,na.rm=T)
  mean_y=mean(y,na.rm=T)
  sd_y=sqrt(var_y)
  
  priorVARe=0.5*var_y
  priorVARg=0.25*var_y
  priorVARb=0.5*sd_y
  priorVARh=0.5*sd_y
  
  default_inits=list(
    list(mu=mean_y, g=rep(0,ng), b=rep(0,ng), h=rep(0,nh), var_e=priorVARe, var_g=priorVARg, var_b=priorVARb, var_h=priorVARh)
  )
  ##inits should be a list of lists, each list cannot have names, 
  ##because rjags will use whether it has names to determin whether it is a list of lists, 
  ##or a single list with variable names. That's bad. 
  #if (!all(c("mu","g","b","h","var_e","var_g","var_h","var_h") %in% names(inits[[i]]))){
  #  stop("mu, g, b, h, var_e, var_g, var_h, var_h for initial values")
  #}
  if(!missing(inits) && !is.null(inits)) {
    if (!is.list(inits)) {
      return("inits parameter must be a list")
    }
    if(length(inits)!=nchain){stop("Number of inital values must be the same as the number of chains")}
    #check names for inits
    inames <- sapply(inits,names)
    if (any(is.null(inames) | nchar(inames) == 0)) {
      return("variable names must be supplied for the initial values")
    }
    null.inits <- sapply(inits, is.null)
    wh.null=which(null.inits)
    n.null=length(wh.null)
  }else{
    inits=default_inits;
    if(nchain>=2)wh.null=c(2:nchain)
    n.null=nchain-1;
  }
  if(n.null>0){
    for(i in 1:n.null){
      set.seed(wh.null[i])
      var_e=runif(1,0.5,2)*priorVARe
      var_g=runif(1,0.5,2)*priorVARg
      var_b=runif(1,0.5,2)*priorVARb
      var_h=runif(1,0.5,2)*priorVARh
      #sample initial values from scaled-inverse-chisquare distribution. But it can sample very rare values.        #df=10
      #var_e=1/rchisq(1,df)*priorVARe*(df+2)	 #same as var_e2=1/rgamma(1,df/2,priorVARe*(df+2)/2)
      #var_g=1/rchisq(1,df)*priorVARg*(df+2)
      #var_b=1/rchisq(1,df)*priorVARb*(df+2)
      #var_h=1/rchisq(1,df)*priorVARh*(df+2)
      inits[[(wh.null[i])]]=list(mu=rnorm(1,mean_y,priorVARe), g=rnorm(ng, 0, priorVARg), b=rnorm(ng, 0, priorVARb), h=rnorm(nh,0,priorVARh), var_e=var_e, var_g=var_g, var_b=var_b, var_h=var_h)
    }            		
  }
  
  return(inits)
}

	
	
	
