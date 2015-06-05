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
  h=h+mu
  IDE=yhat$IDE
  IDL=yhat$IDL
  yhat=yhat$yhat
  y=y[,3]
  uniqIDL=sort(as.numeric(unique(IDL)))
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
