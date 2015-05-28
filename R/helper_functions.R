####################################################################################
### function to plot
####################################################################################
plotYhatvsE=function(yhat=NULL,h=NULL,mu=NULL,g=NULL,b=NULL,IDL=NULL,IDE=NULL,y=NULL,main=NULL){
	#yhat is a dataframe with IDL,IDE and yhat
	# IDL and IDE are the IDLiety index for the predicted yhat values
	#either provide yhat and h or provide h,mu,g,b,IDL,IDE.
	if(is.null(yhat)){
	yhat=data.frame(aggregate(g[IDL]+(1+b[IDL])*h[IDE],by=list(IDL,IDE),mean))
	y=data.frame(aggregate(y,by=list(IDL,IDE),mean))
	colnames(yhat)=c("IDL","IDE","yhat")
    if(!is.null(mu)){yhat$yhat=yhat$yhat+mu}
    }
    h=h+mu
    IDE=yhat$IDE
	IDL=yhat$IDL
	yhat=yhat$yhat
	y=y[,3]
	uniqIDL=unique(IDL)
	n.IDL=length(uniqIDL)
	plot(c(min(c(y,yhat)),max(c(y,yhat)))~ c(min(h),min(h)+(max(h)-min(h))*1.05),type="n",xlab="Environment values",ylab="Variety performance",main=main)
	cols=NULL
	pchs=NULL
	for(i in 1:n.IDL){	
		col=i+1
		pch=1
		cols=c(cols,col)
		pchs=c(pchs,pch)
		IDLi=uniqIDL[i]
		wh.i=which(IDL==IDLi)
		yhat.i=yhat[wh.i]
		y.i=y[wh.i]
		IDE.i=IDE[wh.i]
		p.i=h[IDE.i]
		order.E.i=order(p.i)
		sortp.i=sort(p.i)
		lines(yhat.i[order.E.i]~sortp.i,col=col,type="l")
		points(y.i[order.E.i]~sortp.i,col=col,pch=pch)
		whmax=which.max(p.i)
		#text(x=min(p.i)+1.05*(max(p.i)-min(p.i)),y=y.i[whmax],labels=IDLi,col=col)
	}
	sorth=sort(h[unique(IDE)])
	lines((sorth) ~ sorth, lty=2,col=1)
	legend("bottomright",legend=c(paste("variety",uniqIDL), "slope = 1"),lty=c(rep(1,n.IDL),2),col=c(cols,1))
}

summaryplot=function(IDL,IDE,realizedValue,postMean,LSvalue,plotdir,samps){
	n.IDL=length(unique(IDL))
    n.IDE=length(unique(IDE))
    extend.IDL=rep(1:n.IDL,each=n.IDE)
    extend.IDE=rep(1:n.IDE,n.IDL)
   # predict.group=setdiff(paste(extend.IDL,extend.IDE,sep="_"),paste(IDL,IDE,sep="_"))
   # if(length(predict.group)>0){predict.IDE=as.numeric(gsub("_.*","",predict.group))
   # predict.IDL=as.numeric(gsub(".*_","",predict.group))}

	pdf("plot.pdf")
	plot(realizedValue$h,postMean$h)
    plot(realizedValue$b,postMean$b)
	plot(realizedValue$g,postMean$g)
	plot(realizedValue$h,LSvalue$h)
	plot(realizedValue$b,LSvalue$b)
	plot(realizedValue$g,LSvalue$g)
	plotYhatvsE(h=postMean$h,mu=postMean$mu,g=postMean$g,b=postMean$b,IDL=IDL,IDE=IDE,main="baysian")
	plotYhatvsE(h=LSvalue$h,mu=LSvalue$mu,g=LSvalue$g,b=LSvalue$b,IDL=IDL,IDE=IDE,main="FWlm")
	#plotYhatvsE(LSvalue$yhat,realizedValue$h,main="FWlm against realized value of E")
	#plotYhatvsE(postMean$yhat,realizedValue$h,main="baysian against realized value of E")
	plotYhatvsE(h=realizedValue$h,mu=realizedValue$mu,g=realizedValue$g,b=realizedValue$b,IDL=IDL,IDE=IDE,main="realized variety values against realized E ")
	plotYhatvsE(h=postMean$h,mu=postMean$mu,g=postMean$g,b=postMean$b,IDL=extend.IDL,IDE=extend.IDE,main="baysian extend")
	plotYhatvsE(h=LSvalue$h,mu=LSvalue$mu,g=LSvalue$g,b=LSvalue$b,IDL=extend.IDL,IDE=extend.IDE,main="FWlm extend")
	plotYhatvsE(h=realizedValue$h,mu=realizedValue$mu,g=realizedValue$g,b=realizedValue$b,IDL=extend.IDL,IDE=extend.IDE,main="realized variety values against realized E extend ")
#trace plot
	if(!is.null(samps)){
	require(coda)
	plot(samps)}
	dev.off()
	system(paste("open plot.pdf"))
}

##function to produce summary correlations

corYhat=function(Param1=NULL,Param2=NULL,IDL,IDE){
	
	getyhat=function(Param,IDL,IDE){
	IDL=as.integer(IDL)
	IDE=as.integer(IDE)	
	yhat=data.frame(aggregate(Param$g[IDL]+(1+Param$b[IDL])*Param$h[IDE],by=list(IDL,IDE),mean))[,3]
    if("mu" %in% names(Param)){yhat=yhat+Param$mu}
    return(yhat)
    }
    return(cor(getyhat(Param1,IDL,IDE),getyhat(Param2,IDL,IDE)))

}


summaryCor=function(IDL,IDE,realizedValue,predictedValue){
	corr=rep(0,6)
	n.IDL=length(unique(IDL))
    n.IDE=length(unique(IDE))
    extend.IDL=rep(1:n.IDL,each=n.IDE)
    extend.IDE=rep(1:n.IDE,n.IDL)
    if(length(unique(paste(IDL,IDE)))<length(extend.IDL)){
    	predict.group=setdiff(paste(extend.IDL,extend.IDE,sep="_"),paste(IDL,IDE,sep="_"))
    	predict.IDL=as.integer(gsub("_.*","",predict.group))
   		predict.IDE=as.integer(gsub(".*_","",predict.group))
		corr[5]=corYhat(Param1=realizedValue,Param2=predictedValue,IDL=predict.IDL,IDE=predict.IDE)  
    	corr[6]=corYhat(Param1=realizedValue,Param2=predictedValue,IDL=extend.IDL,IDE=extend.IDE) 
    	}
	
	#dat is a data.frame with the IDL and IDE combinations in the data set
	corr[1]=cor(realizedValue$b,predictedValue$b)
    corr[2]=cor(realizedValue$g,predictedValue$g)
    corr[3]=cor(realizedValue$h,predictedValue$h)
    corr[4]= corYhat(Param1=realizedValue,Param2=predictedValue,IDL=IDL,IDE=IDE)
	names(corr)=c("b","g","h","yhat","predicted yhat only","extended yhat")
	return(corr)	
}




fitmodel=function(y,IDL,IDE,A=NULL,Ainv=NULL,model,nIter,burnIn,thin,seed=NULL,savedir=".",realizedValue=NULL){
		corr=NULL
		for(modeli in model){
			if(modeli=="lm"){
				predictedValue=FWlm(y,IDL,IDE)
				}else
	    		if(modeli=="Gibbs"){
	    			predictedValue=GibbsFW(y=y,IDL=IDL,IDE=IDE,nIter=nIter,burnIn=burnIn,thin=thin,A=A,seed=seed,savedir=savedir)[[1]];
	    			}else
	    			if(modeli=="jags"){
	    				predictedValue = jagsFW(y=y,IDL=IDL,IDE=IDE,Ainv=Ainv,burnIn=burnIn,nIter=nIter,thin=thin,n.adapt=burnIn,seed=seed,savedir=savedir)[[1]];
	    				}else{
	    				error("no model:",modeli)
	    				}
	    	if(!is.null(realizedValue)){			
				corr=cbind(corr,summaryCor(IDL,IDE,realizedValue,predictedValue));
			}
		}
		colnames(corr)=model
		return(corr)
	}




#########################################################
#simulate data
#########################################################
SimuData=function(parameters,savedir,ub="ub_halfVAR",pro.missing=0.5,runModels=T,burnIn=1000,nIter=5000,thin=thin,model){
	#model can be jags, Gibbs, lm
#parameters is a list with var_g,var_b,var_h,var_e,ng,nh,nrep,mu,and(or) A.
	if(!file.exists(savedir))dir.create(savedir,recursive=T)
	for(i in 1:length(parameters)){
		assign(names(parameters)[i],parameters[[i]])
	}
		if(!is.null(parameters$A)){	
		save(A,file=file.path(savedir,"A.rda"))
		g=mvrnorm(n=1,mu=rep(0,ng),Sigma=A*var_g)
		b=mvrnorm(n=1,mu=rep(0,ng),Sigma=A*var_b)
		if("jags"%in%model){Ainv=solve(A);save(Ainv,file=file.path(savedir,"Ainv.rda"))}
	}else{	
		A=NULL
		Ainv=NULL
		g=rnorm(ng,0,sd=sqrt(var_g))
		b=rnorm(ng,0,sd=sqrt(var_b))
	}
	
	h=rnorm(nh,0,sd=sqrt(var_h));
	e=rnorm(ng*nh*nrep,0,sd=sqrt(var_e));	
	#line effect
	IDL=rep(c(1:ng),each=nh*nrep)
	IDE=rep(rep(c(1:nh),each=nrep),ng)
	y=mu+g[IDL]+h[IDE]+b[IDL]*h[IDE]+e

	dat=data.frame(y)
	dat$IDL=IDL
	dat$IDE=IDE
	names(g)=c(1:ng)
	names(h)=c(1:nh)
	names(b)=c(1:ng)
 	realizedValue=list(mu=mu,g=g,h=h,b=b,var_g=var(g),var_p=var(h),var_b=var(b),var_e=var(e))
	save(realizedValue,file=file.path(savedir,"realizedValue.rda"))
	if(!file.exists(file.path(savedir,"balance"))) dir.create(file.path(savedir,"balance"))
	save(dat,file=file.path(savedir,"balance/dat.rda"))
	if(runModels==T) {
		balance.sum=fitmodel(y=y, IDL=IDL, IDE=IDE,A=A,Ainv=Ainv, nIter=nIter, burnIn=burnIn, thin=thin, model=model, seed=NULL,savedir=file.path(savedir,"balance"),realizedValue=realizedValue)
		colnames(balance.sum)=paste("balance",colnames(balance.sum),sep="_")
	     }
	if(! ub %in% c("ub_halfENV","ub_halfVAR","random")){
		stop("ub must be one of the mothods: ub_halfENV, ub_halfVAR, random")
	}	
	if(ub=="ub_halfENV"){
		#unbalanced data: divide ENV into high and low groups, 
		#first half lines in low environments, second five lines in high environments
		lowE=as.integer(names(sort(h)[1:(nh/2)]))
		highE=as.integer(names(sort(h)[(nh/2+1):nh]))

		for(i in 1:nh){
			if(i %in% lowE) dat=dat[-which(dat$IDL%in%c(1:(ng/2)) & dat$IDE==i),]
			else dat=dat[-which(dat$IDL%in%((ng/2+1):ng) & dat$IDE==i),]
		}
	}
    if(ub=="ub_halfVAR"){
    	#unbalanced data: divide VAR into high and low groups,
    	# high half lines in first env group, low half lines in second  env group
		lowG=as.integer(names(sort(g)[1:(ng/2)]))
		highG=as.integer(names(sort(g)[(ng/2+1):ng]))

		for(i in 1:ng){
			if(i %in% lowG) dat=dat[-which(dat$IDE%in%c(1:(nh/2)) & dat$IDL==i),]
			else dat=dat[-which(dat$IDE%in%((nh/2+1):nh) & dat$IDL==i),]
		}
	}
	#randomly missing a porportion pro.missing of genotypes from each enrionment.
	if(ub=="random"){
		ng.remove=round(ng*pro.missing)
		for(i in 1:nh)dat=dat[-which((dat$IDL %in% sample(1:ng,ng.remove)) & dat$IDE==i),]			
	}

	if(!file.exists(file.path(savedir,ub))) dir.create(file.path(savedir,ub))
	save(dat,file=file.path(savedir,ub,"dat.rda"))
	y=dat$y
	IDL=dat$IDL
	IDE=dat$IDE
		
	if(runModels==T) {		ub.sum=fitmodel(y=y,IDL=IDL,IDE=IDE,A=A,Ainv=Ainv,model=model,nIter=nIter,burnIn=burnIn,thin=thin,seed=NULL,savedir=file.path(savedir,ub),realizedValue=realizedValue)
colnames(ub.sum)=paste(ub,colnames(ub.sum),sep="_")		
		corr=cbind(balance.sum,ub.sum)	
		save(corr,file=file.path(savedir,"corr.rda"))}		
		return(corr)
	}

##return significantly different pairs
array_sigdiff=function(dat){
	#three dimentional array data
	#comparisons was done among columns (the third dimension )
	#different parameters within the same comparison was in rows (second dimension)
	#samples was arranged according to the  first dimension of the array 
	nrecords=dim(dat)[1]
	col_names=dimnames(dat)[[3]]
	ncol=length(col_names)
	row_names=dimnames(dat)[[2]]
	nrow=length(row_names)
	#get pairdiff
	#Order columns by the mean
	origroupOrder=col_names
	ToOrder=order(apply(dat,3,mean),decreasing=T)
	dat=dat[,,ToOrder]
	groups=col_names[ToOrder]
	cmb=combn(groups,2)
	npairs=ncol(cmb)
	pairdiff=array(dim=c(nrecords,nrow,npairs))
	for(i in 1:dim(pairdiff)[1]){
		pairdiff[i,,]=apply(cmb,2,function(a){dat[i,,a[1]]-dat[i,,a[2]]})
		
	}
	dimnames(pairdiff)[1:2]=dimnames(dat)[1:2]

#get sigdiff
sigdiff=array(dim=c(nrow,npairs))
for(i in 1:npairs){
	diff=pairdiff[,,i]
	mean_diff=apply(diff,2,mean)
	sd_diff=apply(diff,2,sd)
	p_value=2*pnorm(q=abs(mean_diff),mean=0,sd=sd_diff/sqrt(nrecords),lower.tail=F)
	sigdiff[,i]=p_value<0.05
	}
	sigdiff=apply(sigdiff,1,function(a)cmb[,a,drop=F])
	names(sigdiff)=row_names
	##The order of sigdiff should be from lower to higher, L1, L2, L3: L2 pushed out and increase a letter, then L3 is pushed out and increase a letter.
	#sigdiff=sigdiff[,ncol(sigdiff):1,drop=F]
#bootlab=tstlab(sigdiff,ng=length(groups),groupnames=groups)
#bootlab=bootlab[origroupOrder]
return(sigdiff)
}


