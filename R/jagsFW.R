##this is the baysian implementation with rjags
jagsFW=function(y,IDL,IDE,Ainv=NULL,inits=NULL,nchain=1,burnIn=1000,nIter=5000,thin=1,savedir=".",seed=NULL,n.adapt=0){
n=length(y)
ng=length(unique(IDL))
nh=length(unique(IDE))
data=list(y=y,ng=ng,nh=nh,n=n,IDL=IDL,IDE=IDE,
Vy=var(y))
if(is.null(Ainv)){
modelfile="IaIb.txt"
cat("model
{ 
df<-5
dfg<-5
dfp<-5
dfb<-5
S<-0.5*Vy*(df+2)  
Sg<-0.25*Vy*(dfg+2)
Sb<-0.5*Vy*(dfb+2)  
Sp<-0.5*Vy*(dfp+2)

for (i in 1 : n) {
      y[i] ~ dnorm(mu+g[IDL[i]]+h[IDE[i]]*(1+b[IDL[i]]),tau_e)
      #theta[i]<-h[IDE[i]]*(1+b[IDL[i]])
      }
for ( i in 1:ng){
g[i] ~ dnorm(0,tau_g)
b[i] ~ dnorm(0,tau_b)
   }
   
for (i in 1:nh){
 h[i]~ dnorm(0,tau_h)
}       
##this should be priors, Bugs will use MCMC to sample from the unnormalized posteria.  
   mu ~ dunif(-1E-5,1E05)  
   tau_g ~ dgamma(dfg/2,Sg/2)
   tau_b ~ dgamma(dfb/2,Sb/2)
   tau_h ~ dgamma(dfp/2,Sp/2)
   tau_e ~ dgamma(df/2,S/2)
var_g <- 1/tau_g
var_b<- 1/tau_b
var_h<-1/tau_h
var_e<-1/tau_e
   }",file=modelfile)	
}else{
	data$Ainv=Ainv;
	data$g0=rep(0,ng);
	data$b0=rep(0,ng);
	modelfile="CaCb.txt"
	cat('model
{ 
df<-5
dfg<-5
dfp<-5
dfb<-5
S<-0.5*Vy*(df+2)  
Sg<-0.25*Vy*(dfg+2)
Sb<-0.125*Vy*(dfb+2)  
Sp<-0.125*Vy*(dfp+2)

for (i in 1 : n) {
      y[i] ~ dnorm(mu+g[IDL[i]]+h[IDE[i]]*(1+b[IDL[i]]),tau_e)
      #theta[i]<-h[IDE[i]]*(b[IDL[i]])
      #yhat[i]<-mu+g[IDL[i]]+h[IDE[i]]+theta[i]
      }
      
g[1:ng] ~ dmnorm(g0,Ainv*tau_g)
b[1:ng] ~ dmnorm(b0,Ainv*tau_b)      
   
for (i in 1:nh){
 h[i]~ dnorm(0,tau_h)
}       
##this should be priors, Bugs will use MCMC to sample from the unnormalized posteria.  
   mu ~ dunif(-1E+03,1E+03)  
   tau_g ~ dgamma(dfg/2,Sg/2)
   tau_b ~ dgamma(dfb/2,Sb/2)
   tau_h ~ dgamma(dfp/2,Sp/2)
   tau_e ~ dgamma(df/2,S/2)
   
   var_g <- 1/tau_g
var_b<- 1/tau_b
var_h<-1/tau_h
var_e<-1/tau_e
   }
',file=modelfile)
	}
#save(data,file="data.rda");setwd(datadir)

parameters<-c("mu","g","b","h","var_g","var_b","var_h","var_e")

############################################# 
# initialize
########################################################################################## 
inits=initialize(y,ng=ng,nh=nh,model="jags",inits=inits,jags.seed=seed,nchain=nchain)


jags.m<-jags.model(file=modelfile,data=data,inits=inits,n.chains=length(inits),n.adapt=n.adapt)
#list.samplers(jags.m)
samps<-coda.samples(jags.m,parameters,n.iter=nIter,thin=thin)
save(samps,file=file.path(savedir,paste("samps.rda",sep="")))
ans=list()
for(i in 1:nchain){
	
	xi=samps[[i]][round(burnIn/thin):round(nIter/thin),,drop=F]
	tmp=apply(xi,2,mean)
	postMean=list()
	postMean$h=tmp[paste("h[",1:data$nh,"]",sep="")]
	postMean$b=tmp[paste("b[",1:data$ng,"]",sep="")]
	postMean$g=tmp[paste("g[",1:data$ng,"]",sep="")]
	postMean$mu=tmp['mu']
	#save(postMean,file=file.path(savedir,paste("postMeanInit",i,".rda",sep="")))
ans[[i]]=postMean
}
names(ans)=paste("Init",c(1:length(inits)),sep="")
file.remove(modelfile)
save(ans,file=file.path(savedir,"postMean_jags.rda"))
return(ans)
#if(length(inits)>1){return(gelman.diag(samps))}
}






