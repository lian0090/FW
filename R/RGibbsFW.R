##This is the Gibbs sampler completely written in R, it is slow!
RGibbsFW=function(y,IDL,IDE,Ainv=NULL,nIter=3000,burnIn=1000,thin=5,savedir=".",seed=NULL){
if(!is.null(seed)) set.seed(seed)	
library(MASS)
current.dir=getwd()	
if(!file.exists(savedir)){dir.create(savedir)}	
setwd(savedir)

#hyper parameters:
df=dfg=5
dfh=dfb=5
S<-0.5*var(y)*(df+2)  #S is the scale times df
Sg<-0.25*var(y)*(dfg+2)
Sb<-0.5*var(y)*(dfb+2)  #this affects the sampler
Sh<-0.5*var(y)*(dfh+2)  #this affects the sampler

############################################# 
# IDs and vector sizes
#############################################
ng=length(unique(IDL))
nSamples<-nIter-burnIn+1
n=length(y)
nh=length(unique(IDE))

############################################# 
# initialize
########################################################################################## 

 
inits1=initialize(y,ng=ng,nh=nh,model="Gibbs",inits=NULL,seed=seed)[[1]]
for(k in 1:length(inits1)){
	assign(names(inits1)[k],inits1[[k]])

}
e<-y-mu-g[IDL]-h[IDE]-b[IDL]*h[IDE]


############################################# 
# Posteria storage
#############################################
post_mu<-0
post_g=rep(0,ng)
post_b=rep(0,ng)
post_h=rep(0,nh)
post_var_e<-0
post_var_g=0
post_var_b=0;
post_var_h=0


############################################# 
# Sampling
#############################################
for(i in 1:nIter){
#cat(i,"var_g=",var_g, "var_h=",var_h,"var_b=",var_b,"h[1]=",h[1],"b[1]=",b[1],mu,"\n")

#environment effect
X<-1+b[IDL]
e<-e+X*h[IDE]
rhs<-tapply(INDEX=IDE,FUN=sum,X=(X*e))/var_e 
C<-tapply(INDEX=IDE,FUN=sum,X=(X*X))/var_e  + 1/var_h  
h<-rnorm(n=nh,sd=sqrt(1/C),mean=rhs/C)
X<-1+b[IDL]
e<-e-X*h[IDE]

#effect of b 
X<-h[IDE]
e=e+X*b[IDL]
if(is.null(Ainv)){
	rhs<-tapply(INDEX=IDL,FUN=sum,X=(X*e))/var_e    
   C<-tapply(INDEX=IDL,FUN=sum,X=(X*X))/var_e  + 1/var_b  
   b<-rnorm(n=ng,sd=sqrt(1/C),mean=rhs/C)
  }else{
  rhs<-tapply(INDEX=IDL,FUN=sum,X=(X*e))/var_e    
  C<-diag(tapply(INDEX=IDL,FUN=sum,X=(X*X)))/var_e  + Ainv/var_b
  Cinv=solve(C)
  b<-mvrnorm(n=1,mu=Cinv%*%rhs,Sigma=Cinv)
 } 
e=e-b[IDL]*h[IDE]
#sample g
e=e+g[IDL]
if(is.null(Ainv)){
Vghat=(aggregate(x=rep(1,n),FUN=sum,by=list(IDL))[,2]+var_e/var_g)^{-1}*var_e
ghat=Vghat/var_e*aggregate(x=e,FUN=sum,by=list(IDL))[,2]
g<-rnorm(n=ng,mean=ghat,sd=sqrt(Vghat))
}else
{
C=(diag(aggregate(x=rep(1,n),FUN=sum,by=list(IDL))[,2])/var_e+Ainv/var_g)
rhs=aggregate(x=e,FUN=sum,by=list(IDL))[,2]/var_e
Cinv=solve(C)
g<-mvrnorm(n=1,mu=Cinv%*%rhs,Sigma=Cinv)
}
e=e-g[IDL]


#sample one g and b together

# for(k in 1:ng){
    # which.k=which(IDL==k)
    # e.k=e[which.k]
    # #sample b
	# h.k=(h[IDE])[which.k]
	# e.k=e.k+h.k*b[k]
   # C=1/var_b+sum(h.k^2)/var_e
    # rhs=sum(h.k*e.k)/var_e
    # b[k]=rnorm(n=1,mean=rhs/C,sd=sqrt(1/C))
   # e.k=e.k-h.k*b[k]
    
	# #sample g
	# e.k=e.k+g[k]
    # C=1/var_g+length(which.k)/var_e
   # rhs=sum(e.k)/var_e
    # g[k]=rnorm(n=1,mean=rhs/C,sd=sqrt(1/C))
    # e.k=e.k-g[k]
	
	# e[which.k]=e.k

	    # }



#var_e
SS<-sum(e^2)+S
DF<-n+df
var_e<-SS/rchisq(df=DF,n=1)

#var_h
SS<-sum(h^2)+Sh
DF<-nh+dfh
var_h<-SS/rchisq(df=DF,n=1)


#var_g
SS<-ifelse(is.null(Ainv),sum(g^2)+Sg,as.numeric(t(g)%*%Ainv%*%g)+Sg)
DF<-ng+dfg
var_g<-SS/rchisq(df=DF,n=1)

#var_b
SS<-ifelse(is.null(Ainv),sum(b^2)+Sb,as.numeric(t(b)%*%Ainv%*%b)+Sb)
DF<-ng+dfb
var_b<-SS/rchisq(df=DF,n=1)


# intercept
#mu=realizedValue$mu
e<-e+mu
mu<-rnorm(1,sd=sqrt(var_e/n),mean=mean(e))
e=e-mu #e need to be updated here because this is needed for updating the corrected y.



if(i%%thin==0){
	append=ifelse(i>thin,T,F)
	write(mu,file="mu.dat",append=append)
	write(b,file="b.txt",ncolumns=ng,append=append,sep=",")
	write(g,file="g.txt",ncolumns=ng,append=append,sep=',')
	write(h,file="p.txt",ncolumns=nh,append=append,sep=",")
	write(var_g,file="var_g.dat",append=append)
	write(var_h,file="var_h.dat",append=append)
	write(var_b,file="var_b.dat",append=append)
	write(var_e,file="var_e.dat",append=append)
}

if(i>=burnIn){
	post_mu=post_mu+mu/nSamples
	post_g=post_g+g/nSamples
	post_h=post_h+h/nSamples
	post_b=post_b+b/nSamples
	post_var_g=post_var_g+var_g/nSamples
	post_var_h=post_var_h+var_h/nSamples
	post_var_e=post_var_e+var_e/nSamples
	post_var_b=post_var_b+var_b/nSamples
}
}
postMean=list(post_mu,post_g,post_h,post_b,post_var_g,post_var_h,post_var_b,post_var_e)
names(postMean)=c("mu","g","h","b","var_g","var_h","var_b","var_e")
postMean=list(Init1=postMean)
#save(postMean,file="postMean.rda")
setwd(current.dir)
return(postMean)
}