lmFW=function(y,VAR,ENV,VARlevels=NULL,ENVlevels=NULL,savedir="."){

   if(!file.exists(savedir)){dir.create(savedir)}	


  IDEL=getIDEL(VAR,ENV,VARlevels,ENVlevels)
  for(i in 1:length(IDEL)){
  	assign(names(IDEL)[i],IDEL[[i]])
  }


  IDE=IDEL$IDE

  IDL=IDEL$IDL

  VARlevels=IDEL$VARlevels

  ENVlevels=IDEL$ENVlevels

  fVAR=IDEL$fVAR

  h=tapply(y,INDEX=IDE,mean)-mean(y)

  n.var=length(VARlevels)

  ZX=sweep(model.matrix(~fVAR-1),1,h[IDE],"*")

  lm1=lm(y~-1+fVAR+ZX)

  gb=coef(lm1)

  g=gb[paste("fVAR",VARlevels,sep="")]

  b=gb[paste("ZXfVAR",VARlevels,sep="")]-1

  names(g)=VARlevels

  names(b)=VARlevels

  names(h)=ENVlevels

  LSvalue=setFW(g=g,b=b,h=h,y=y,VAR=VAR,ENV=ENV,IDL=IDL,IDE=IDE,VARlevels=VARlevels,ENVlevels=ENVlevels,mu=0) 

  class(LSvalue)=c("FW","list")

  save(LSvalue,file=file.path(savedir,"LSvalue.rda"))

  return(LSvalue)

}




	
