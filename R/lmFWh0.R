lmFWh0=function(y,VAR,ENV){
  #if genotype or environment is completely missing for a GxE combination, the predicted value of  y is still NA.
  VAR = factor(VAR)
  ENV = factor(ENV)
  IDL = as.numeric(VAR)
  IDE = as.numeric(ENV)
  VARlevels = levels(VAR)
  ENVlevels = levels(ENV)
  n.var=length(VARlevels)
  n.env=length(ENVlevels)
  ##step 1 obtain environment effect (with sum contrast)
 # lm0=lm(y~-1+fENV+fVAR)
 # h=coef(lm0)[paste("fENV",ENVlevels,sep="")]
#  h=h-mean(h,na.rm=T)
  fVARc=VAR; attr(fVARc, "contrasts") <- contr.sum(n.var) 
  fENVc=ENV; attr(fENVc,"contrasts")<-contr.sum(n.env)
  mf=model.frame(y~fENVc+fVARc)
  lm0=lm(mf)
  h=coef(lm0)[c(2:n.env)]
  h=c(h,-sum(h,na.rm=T))
 # h=tapply(y,INDEX=IDE,function(a)mean(a,na.rm=T))-mean(y,na.rm=T) 
  g=rep(0,n.var)
  b=rep(0,n.var)
  var_e=rep(0,n.var)
  df=rep(0,n.var)
  for(i in 1:n.var){
   whVar=which(IDL==i)
   lmi=lm(y[whVar]~h[IDE[whVar]])
   sum_lmi=summary(lmi)
   g[i]=coef(lmi)[1]
   b[i]=coef(lmi)[2]-1
   df[i]=sum_lmi$df[2]
   var_e[i]=(sum_lmi$sigma)^2
  }
 g=matrix(g)
 b=matrix(b)  
 h=matrix(h)
 rownames(g)=VARlevels
 rownames(b)=VARlevels
 rownames(h)=ENVlevels
 yhat=matrix(g[IDL,]+(1+b[IDL,])*h[IDE,]) #do not use VAR or ENV for index, because they might be factor
 ##this var_e will be exactly the same as if fitting all observations in a single linear model
 var_e_weighted=sum(var_e*df)/sum(df)
 
 LSvalue=list(y=y, whichNa = which(is.na(y)), VAR = VAR, ENV = ENV, VARlevels=VARlevels,ENVlevels=ENVlevels, mu = 0, g = g, b = b, h = h, yhat = yhat,var_e=var_e, var_e_weighted=var_e_weighted) 
 class(LSvalue)=c("FW","list")
 return(LSvalue)
}





