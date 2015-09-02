lmFW=function(y,VAR,ENV,VARlevels=NULL,ENVlevels=NULL){
  #if genotype or environment is completely missing for a GxE combination, the predicted value of  y is still NA.
  IDEL=getIDEL(VAR,ENV,VARlevels,ENVlevels)
  for(i in 1:length(IDEL)){
    assign(names(IDEL)[i],IDEL[[i]])
  } 
  h=tapply(y,INDEX=IDE,function(a)mean(a,na.rm=T))
  n.var=length(VARlevels)
  g=rep(0,n.var)
  b=rep(0,n.var)
  for(i in 1:n.var){
  	whVar=which(IDL==i)
  	lmi=lm(y[whVar]~h[IDE[whVar]])
  	g[i]=coef(lmi)[1]
  	b[i]=coef(lmi)[2]-1
  }
  g=matrix(g)
  b=matrix(b)  
  h=matrix(h)
  rownames(g)=VARlevels
  rownames(b)=VARlevels
  rownames(h)=ENVlevels
  yhat=matrix(g[VAR,]+(1+b[VAR,])*h[ENV,])
  LSvalue=list(y=y, whichNa = which(is.na(y)), VAR = VAR, ENV = ENV, VARlevels = VARlevels, ENVlevels = ENVlevels, g = g, b = b, h = h, yhat = yhat) 
  class(LSvalue)=c("FW","list")
  return(LSvalue)
}





