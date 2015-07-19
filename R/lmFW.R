lmFW=function(y,VAR,ENV,VARlevels=NULL,ENVlevels=NULL){

  #if genotype or environment is completely missing for a GxE combination, the predicted value of  y is still NA.
  IDEL=getIDEL(VAR,ENV,VARlevels,ENVlevels)
  for(i in 1:length(IDEL)){
  	assign(names(IDEL)[i],IDEL[[i]])
  }

    
  h=tapply(y,INDEX=IDE,function(a)mean(a,na.rm=T))-mean(y,na.rm=T)

  n.var=length(VARlevels)

  ZX=sweep(model.matrix(~fVAR-1),1,h[IDE],"*")

  lm1=lm(y~-1+fVAR+ZX)

  gb=coef(lm1)

  g=matrix(gb[paste("fVAR",VARlevels,sep="")])

  b=matrix(gb[paste("ZXfVAR",VARlevels,sep="")]-1)

  h=matrix(h)
  
  rownames(g)=VARlevels

  rownames(b)=VARlevels

  rownames(h)=ENVlevels

  yhat=matrix(g[VAR,]+(1+b[VAR,])*h[ENV,])

 
  LSvalue=list(y=y,whichNa=which(is.na(y)),VAR=VAR,ENV=ENV,VARlevels=VARlevels,ENVlevels=ENVlevels,mu=0,g=g,b=b,h=h,yhat=yhat) 

  class(LSvalue)=c("FW","list")


  return(LSvalue)

}




	
