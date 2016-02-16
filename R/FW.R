#a wrapper for GibbsFW and lmFW
#VAR and ENV coersed to factors
FW = function(y, VAR, ENV,  method = c("OLS", "Gibbs")[2], A = NULL, H = NULL, saveAt = "", nIter = 5000, 
	burnIn = 3000, thin = 5, dfe = 5, dfg = 5, dfh = 5, dfb = 5, priorVARe = NULL, priorVARg = NULL, priorVARb = NULL, priorVARh = NULL, nchain = 1, 
	seed = NULL, inits = NULL, saveVAR = c(1:2), saveENV = c(1:2)) {
	model = "h0"

	if (nIter <= burnIn) 
		stop("nIter must be larger than burnIn\n")

	if(any(is.na(VAR)))stop("NA  in VAR is not allowed\n")
  
  if(any(is.na(ENV)))stop("NA in ENV is ont allowed\n")
	#coerse to factors
   
	if (method == "OLS") {
		predictedValue = lmFWh0(y, VAR, ENV)
	}
	if (method == "Gibbs") {
		predictedValue = GibbsFWh0(y = y, VAR = VAR, ENV = ENV,  nIter = nIter, burnIn = burnIn, 
			thin = thin, dfe = dfe, dfg = dfg, dfh = dfh, dfb = dfb, priorVARe = priorVARe, priorVARg = priorVARg, priorVARb = priorVARb, priorVARh = priorVARh, 
			A = A, H = H, nchain = nchain, seed = seed, inits = inits, saveAt = saveAt,saveVAR=saveVAR,saveENV=saveENV)
	}


	return(predictedValue)

}


####################################################################################
### function to plot FW object
####################################################################################
plot.FW = function(x, plotVAR=NULL, chain=1, ENVlabel="split", ...) {
  ##first argument name must be x to be consistent with S3 generic methods
  if(is.null(chain))chain=1
 y=x$y
 VARlevels=x$VARlevels
 ENVlevels=x$ENVlevels
 VAR=factor(x$VAR,levels=VARlevels)
 ENV=factor(x$ENV,levels=ENVlevels)
 IDL=as.numeric(VAR)
 IDE=as.numeric(ENV)
 yhat=x$yhat[,chain]
 mu=x$mu[chain]
 g=x$g[,chain]
 names(g)=VARlevels
 b=x$b[,chain]
 names(b)=VARlevels
 h=x$h[,chain]
 names(h)=ENVlevels
  
  if (!is.null(plotVAR)) {
    if (is.integer(plotVAR)) {
      plotIDL = plotVAR
    } else if (is.character(plotVAR)) {
      plotIDL = match(plotVAR,VARlevels,nomatch=0) 
    }
  }else{
    plotIDL=c(1:length(VARlevels))
  } 
  n.plotVAR = length(plotIDL)
  oripar = par()$mar
  if(ENVlabel==FALSE){
  par(xpd = T, mar = oripar + c(1.5, 0, 0, 5))
  }else{
  par(xpd = T, mar = oripar + c(0, 0, 0, 5))  
  }
  range.h = range(h, na.rm = T)
 #set range.y with ys in plotIDL
 y.plot=y[(IDL%in%plotIDL)]
 range.y = range(y.plot, na.rm = T)
  
  #default valus for ...
  args1=list(xlab="Environment effect",ylab="Variety performance",bty="l",type="n")
  inargs<-list(...)
  args1[names(inargs)]<-inargs
  
  if(ENVlabel==FALSE){
  do.call(plot,c(list(formula=range.y ~ range.h),args1))
  }else{
    args2=args1
    args2$xlab=""
    do.call(plot,c(list(formula=range.y ~ range.h), args2))
    args2=args1
    args2$main=""
    args2$ylab=""
    args2$line=4
    args2$type=NULL
    do.call(title,args2)    
  }
  sorth = sort(h)
  
  sorth1 = sorth[seq(1, length(h), by = 2)]
  sorth2 = sorth[seq(2, length(h), by = 2)]
  if(ENVlabel!=FALSE){
    if(ENVlabel=="split"|ENVlabel==TRUE){
      axis(side = 1, at = sorth1, labels = names(sorth1), line = 2)
      axis(side = 3, at = sorth2, labels = names(sorth2), line = 1)
    }else if(ENVlabel=="bottom"){
      axis(side = 1, at=sorth,labels=names(sorth),line=2)
    }else if(ENVlabel=="top"){
      axis(side=3,at=sorth,labels=names(sorth),line=1)
    }else{
      cat("ENVlabel must be TRUE, FALSE, \"split\",\"top\" or \"bottom\" \n")
    }
  }
  cols = NULL
  pchs = NULL
 
  usr <- par("usr")
  col=1
  for (i in plotIDL) {
    col = col + 1
    pch = 1
    cols = c(cols, col)
    pchs = c(pchs, pch)
    IDLi = which(IDL == i)
    y.i = y[IDLi]
    ENV.i=ENV[IDLi]
    #aggregate will not change the factor levels  
    cellMeans=aggregate(y.i,by=list(ENV.i),mean)
    cellMeansIDE=as.numeric(cellMeans[,1])
    cellMeans=cellMeans[,2]
    hi=h[cellMeansIDE]
    clip(range.h[1], range.h[2], range.y[1], range.y[2])
    abline(a=mu+g[i],b=b[i]+1, col = col, ...)
    do.call("clip", as.list(usr))
    points(cellMeans ~ hi, col = col, pch = pch,xpd=T,...)    
  }
  #the intercept denpends on g
  clip(range.h[1], range.h[2], range.y[1], range.y[2])
  abline(a = mean(y.plot, na.rm = T), b = 1, lty = 2, col = 1)
  clip(usr[1], usr[2] + (usr[2]-usr[1])* 10, usr[3], usr[4])
  legend(x = usr[2], y = usr[4], legend = c(VARlevels[plotIDL], "slope = 1"), lty = c(rep(1, n.plotVAR), 2), col = c(cols, 1), bg = "transparent", bty = "n")
  par(mar = oripar)
  do.call("clip", as.list(usr))
}
####################################################################################
### print postMean and FW object
####################################################################################

#print.FW=function(FWobj){
#  cat("FW object\n")
# print(FWobj[c("mu","g","b","h")])
#}


getIDEL = function(VAR, ENV) {
  #coerse to factors
	out$IDL = as.numeric(VAR)
	out$IDE = as.numeric(ENV)
  out$VARlevels = levels(VAR)
	out$ENVlevels = levels(ENV)

	return(out)
}



