#initialize
initialize=function(y,ng,nh,model,inits=NULL,seed=NULL,nchain=1){
	#seed is to set the random seed for Gibbs Sampler for jags. Not for the random seed of inital values.
	var_y=var(y,na.rm=T)
	sd_y=sqrt(var_y)
	if(model=="Gibbs"){
		default_inits=list(
			list(mu=mean(y), g=rep(0,ng), b=rep(0,ng), h=rep(0,nh), var_e=sd_y, var_g=sd_y, var_b=sd_y, var_h=sd_y)
		)
		}
	else if(model=="jags"){
		default_inits=list(
			list(mu=mean(y),g=rep(0,ng),b=rep(0,ng),h=rep(0,nh),tau_e=1/sd_y,tau_g=1/sd_y, tau_b=1/sd_y, tau_h=1/sd_y, .RNG.seed=seed, .RNG.name="base::Mersenne-Twister")
			)
	}else{
		stop("model must either be Gibbs or jags")
		}	
		
       ##inits should be a list of lists, each list cannot have names, 
       ##because rjags will use whether it has names to determin whether it is a list of lists, 
       ##or a single list with variable names. That's bad. 
	          
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
            		if(model=="Gibbs"){
            			inits[[(wh.null[i])]]=list(mu=rnorm(1,mean(y),var_y/2), g=rnorm(ng, 0, var_y/4), b=rnorm(ng, 0, var_y/2), h=rnorm(nh,0,var_y/2), var_e=var_y/2, var_g=var_y/4, var_b=var_y/2, var_h=var_y/2)}else if(model=="jags"){
            			inits[[(wh.null[i])]]=list(mu=rnorm(1,mean(y),var_y/2), g=rnorm(ng, 0, var_y/4), b=rnorm(ng, 0, var_y/2), h=rnorm(nh,0,var_y/2), tau_e=2/var_y, tau_g=4/var_y, tau_b=2/var_y, tau_h=2/var_y, .RNG.seed=seed, .RNG.name="base::Mersenne-Twister")
            	}
            }	
    	}
      
  	return(inits)
 }

	
	
	
