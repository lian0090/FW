#a wrapper for GibbsFW and lmFW
FW=function(y,VAR,ENV,method=c("OLS","Gibbs")[2], saveAt="",nIter=5000,burnIn=3000,thin=5,df=5,dfg=5,dfh=5,dfb=5,priorVARe=NULL,priorVARg=NULL,priorVARb=NULL,priorVARh=NULL,A=NULL,nchain=1,seed=NULL,saveVAR=c(1:2),saveENV=c(1:2)){


if(saveAt==""){
	saveAt=paste(getwd(),"/",sep="");
}

whNA=which(is.na(y))
#if genotype or environment is completely missing for a GxE combination, the predicted value of  y is still NA.

if(method=="OLS"){
	 predictedValue=lmFW(y,VAR,ENV)
}

if(method=="Gibbs"){
	 	 predictedValue=GibbsFW(y=y,VAR=VAR,ENV=ENV,nIter=nIter,burnIn=burnIn,thin=thin,df=df,dfg=dfg,dfh=dfh,dfb=dfb,priorVARe=priorVARe,priorVARg=priorVARg,priorVARb=priorVARb,priorVARh=priorVARh, A=A,nchain=nchain,seed=seed,saveAt=saveAt);	
	
}

return(predictedValue)

}


####################################################################################
### function to plot FW object
####################################################################################
plot.FW=function(FWobj,plotVAR=NULL,main=NULL,chain=1){
  #plotVAR: the variety names for which the plot should be generated
  #default to print the frist chain
  namesFWobj=names(FWobj)
  for(i in 1:length(namesFWobj)){
    assign(namesFWobj[i],FWobj[[i]])
  }

  if(!is.null(plotVAR)){
    if(is.numeric(plotVAR)){
    	plotIDL=plotVAR
    	    	}else if (is.character(plotVAR)){
    	    	plotIDL=which(VARlevels %in% plotVAR)
    	    		
    	    	}
    	
    VARlevels=VARlevels[plotIDL]
    whVAR=which(VAR %in% VARlevels)
    y=y[whVAR]
    VAR=VAR[whVAR]
    ENV=ENV[whVAR]
    yhat=yhat[whVAR,chain]
  }else{
  	yhat=yhat[,chain]
  }
  yhat=aggregate(yhat,by=list(VAR,ENV),function(a)mean(a,na.rm=T))
  y=data.frame(aggregate(y,by=list(VAR,ENV),function(a)mean(a,na.rm=T)))
  colnames(yhat)=c("VAR","ENV","yhat")
  yhat$yhat=yhat$yhat
  VAR=yhat$VAR
  ENV=yhat$ENV
  yhat=yhat$yhat
  y=y[,3]
  n.VAR=length(VARlevels)
  plot(c(min(c(y,yhat),na.rm=T),max(c(y,yhat),na.rm=T))~ c(min(h,na.rm=T),min(h,na.rm=T)+(max(h,na.rm=T)-min(h,na.rm=T))*1.05),type="n",xlab="",ylab="Variety performance",main=main)
 
  h=h[order(h[,chain]),chain]
  
  sorth1=h[seq(1,length(h),by=2)]
  sorth2=h[seq(2,length(h),by=2)]
  axis(side=1,at=sorth1,labels=names(sorth1),line=1,las=2)
  axis(side=3,at=sorth2,labels=names(sorth2),line=1,las=2)
  cols=NULL
  pchs=NULL
  for(i in 1:n.VAR){  
    col=i+1
    pch=1
    cols=c(cols,col)
    pchs=c(pchs,pch)
    VARi=VARlevels[i]
    wh.i=which(VAR==VARi)
    yhat.i=yhat[wh.i]
    y.i=y[wh.i]
    ENV.i=ENV[wh.i]
    p.i=h[ENV.i]
    order.E.i=order(p.i)
    sortp.i=sort(p.i)
    lines(yhat.i[order.E.i]~sortp.i,col=col,type="l")
    points(y.i[order.E.i]~sortp.i,col=col,pch=pch)
    whmax=which.max(p.i)
    #text(x=min(p.i)+1.05*(max(p.i)-min(p.i)),y=y.i[whmax],labels=IDLi,col=col)
    
  }
  abline(a=mean(g[,chain]+mu[chain]),b=1,lty=2,col=1)
  legend("bottomright",legend=c(VARlevels, "slope = 1"),lty=c(rep(1,n.VAR),2),col=c(cols,1))
}
####################################################################################
### print postMean and FW object
####################################################################################

#print.FW=function(FWobj){
#  cat("FW object\n")
 # print(FWobj[c("mu","g","b","h")])
#}


getIDEL=function(VAR,ENV,VARlevels=NULL,ENVlevels=NULL){
  VAR=as.character(VAR)
  ENV=as.character(ENV)
  if(is.null(VARlevels)){
  VARlevels=sort(unique(VAR))
  }
  if(is.null(ENVlevels)){
    ENVlevels=sort(unique(ENV))
  }
  fVAR=factor(VAR,levels=VARlevels,ordered=T)
  fENV=factor(ENV,levels=ENVlevels,ordered=T)
  IDL=as.numeric(fVAR)
  IDE=as.numeric(fENV)
  out=list()
  out$IDL=IDL
  out$IDE=IDE
  out$VARlevels=VARlevels
  out$ENVlevels=ENVlevels
  out$fVAR=fVAR
  out$fENV=fENV
  return(out)
}



