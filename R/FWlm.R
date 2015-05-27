#this is the code for ordinary linear regression
FWlm=function(y,IDL,IDE,savedir="."){
	if(!is.numeric(IDL)){stop("IDL must be numeric index for varieties")}
	if(!is.numeric(IDE)){stop("IDL must be numeric index for envrionments")}
	
	h=tapply(y,INDEX=IDE,mean)-mean(y)
	n.var=length(unique(IDL))
	fVAR=factor(IDL,levels=1:n.var,ordered=T)
	biEj=diag(h[IDE])%*%model.matrix(~fVAR-1)
	lm1=lm(y~-1+fVAR+biEj)
	ab=coef(lm1)
	g=ab[paste("fVAR",c(1:n.var),sep="")]
	b=ab[paste("biEjfVAR",c(1:n.var),sep="")]-1
	LSvalue=list(g=g,b=b,h=h)
	save(LSvalue,file=file.path(savedir,"LSvalue.rda"))
	return(LSvalue)
	
}


