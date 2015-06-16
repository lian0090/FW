####################################################################################
### function to plot FW object
####################################################################################
plot.FW=function(FWobj,plotVAR=NULL,main=NULL){
  #plotVAR: the variety names for which the plot should be generated
  namesFWobj=names(FWobj)
  for(i in 1:length(namesFWobj)){
    assign(namesFWobj[i],FWobj[[i]])
  }
  if(!is.null(plotVAR)){
    plotIDL=which(VARlevels %in% plotVAR)
    VARlevels=VARlevels[plotIDL]
    whIDL=which(IDL %in% plotIDL)
    IDL=IDL[whIDL]
    IDE=IDE[whIDL]
    fitted.values=fitted.values[whIDL]
    y=y[whIDL]
  }
  yhat=aggregate(fitted.values,by=list(IDL,IDE),mean)
  y=data.frame(aggregate(y,by=list(IDL,IDE),mean))
  colnames(yhat)=c("IDL","IDE","yhat")
  yhat$yhat=yhat$yhat+mu
  IDE=yhat$IDE
  IDL=yhat$IDL
  yhat=yhat$yhat
  y=y[,3]
  uniqIDL=sort(as.numeric(unique(IDL)))
  n.IDL=length(uniqIDL)
  plot(c(min(c(y,yhat)),max(c(y,yhat)))~ c(min(h),min(h)+(max(h)-min(h))*1.05),type="n",xlab="",ylab="Variety performance",main=main)
  sorth=sort(h)
  sorth1=sorth[seq(1,length(h),by=2)]
  sorth2=sorth[seq(2,length(h),by=2)]
  axis(side=1,at=sorth1,labels=names(sorth1),line=1,las=2)
  axis(side=3,at=sorth2,labels=names(sorth2),line=1,las=2)
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
  sorth=sort(h)
  abline(a=mean(g+mu),b=1,lty=2,col=1)
  legend("bottomright",legend=c(VARlevels, "slope = 1"),lty=c(rep(1,n.IDL),2),col=c(cols,1))
}
####################################################################################
### print postMean and FW object
####################################################################################
print.postMean=function(postMean){
  for(i in 1:length(postMean)){
    print((postMean[[i]])[c("var_e","var_g","var_b","var_h","mu","g","b","h")] )
  }  
}
print.FW=function(FWobj){
  cat("FW object\n")
  print(FWobj[c("mu","g","b","h")])
}


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


setFW=function(g,b,h,y,VAR,ENV,...){
  eclipse=list(...)
  fitted.values=g[VAR]+(1+b[VAR])*h[ENV]
  ENVmean=aggregate(y,by=list(ENV),mean)
  nameENVmean=ENVmean[,1]
  ENVmean=ENVmean[,2]
  names(ENVmean)=nameENVmean
  ENVmean=ENVmean[names(h)]
  corENVmean=cor(ENVmean,h)
  corfitted=cor(y,fitted.values)
  cor_ymean=get_cor_ymean(g=g,b=b,h=h,y=y,VAR=VAR,ENV=ENV)
  out=list(g=g,b=b,h=h,y=y,VAR=VAR,ENV=ENV,ENVmean=ENVmean,corENVmean=corENVmean,corfitted=corfitted,cor_ymean=cor_ymean,fitted.values=fitted.values,mu=0)
  out=c(out,eclipse)	
  class(out)=c("FW","list")
  
  return(out)
  
}
