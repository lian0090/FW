#a wrapper for GibbsFW and lmFW
FW = function(y, VAR, ENV, VARlevels = NULL, ENVlevels = NULL, method = c("OLS", "Gibbs")[2], A = NULL, H = NULL, saveAt = "", nIter = 5000, 
	burnIn = 3000, thin = 5, dfe = 5, dfg = 5, dfh = 5, dfb = 5, priorVARe = NULL, priorVARg = NULL, priorVARb = NULL, priorVARh = NULL, nchain = 1, 
	seed = NULL, inits = NULL, saveVAR = c(1:2), saveENV = c(1:2)) {
	model = "h0"

	if (nIter <= burnIn) 
		stop("nIter must be larger than burnIn\n")

	if(any(is.na(VAR)))stop("NA  in VAR is not allowed\n")
  
  if(any(is.na(ENV)))stop("NA in ENV is ont allowed\n")

	if (method == "OLS") {
		predictedValue = lmFWh0(y, VAR, ENV, VARlevels = VARlevels, ENVlevels = ENVlevels)
	}
	if (method == "Gibbs") {
		predictedValue = GibbsFWh0(y = y, VAR = VAR, ENV = ENV, VARlevels = VARlevels, ENVlevels = ENVlevels, nIter = nIter, burnIn = burnIn, 
			thin = thin, dfe = dfe, dfg = dfg, dfh = dfh, dfb = dfb, priorVARe = priorVARe, priorVARg = priorVARg, priorVARb = priorVARb, priorVARh = priorVARh, 
			A = A, H = H, nchain = nchain, seed = seed, inits = inits, saveAt = saveAt,saveVAR=saveVAR,saveENV=saveENV)
	}


	return(predictedValue)

}


####################################################################################
### function to plot FW object
####################################################################################
plot.FW = function(x, plotVAR=NULL, chain=1, ENVlabel=T, ...) {
  ##first argument name must be x to be consistent with S3 generic methods
  namesx = names(x)
  for (i in 1:length(namesx)) {
    assign(namesx[i], x[[i]])
  }
  if(is.null(chain))chain=1
  if (!is.null(plotVAR)) {
    if (is.integer(plotVAR)) {
      plotIDL = plotVAR
    } else if (is.character(plotVAR)) {
      plotIDL = which(VARlevels %in% plotVAR)
      
    }
    
    VARlevels = VARlevels[plotIDL]
    whVAR = which(VAR %in% VARlevels)
    y = y[whVAR]
    VAR = VAR[whVAR]
    ENV = ENV[whVAR]
    yhat = yhat[whVAR, chain]
  } else {
    yhat = yhat[, chain]
  }
  yhat = aggregate(yhat, by = list(VAR, ENV), function(a) mean(a, na.rm = T))
  y = data.frame(aggregate(y, by = list(VAR, ENV), function(a) mean(a, na.rm = T)))
  colnames(yhat) = c("VAR", "ENV", "yhat")
  yhat$yhat = yhat$yhat
  VAR = yhat$VAR
  ENV = yhat$ENV
  yhat = yhat$yhat
  y = y[, 3]
  n.VAR = length(VARlevels)
  oripar = par()$mar
  if(!ENVlabel){
  par(xpd = T, mar = oripar + c(1.5, 0, 0, 5))
  }else{
  par(xpd = T, mar = oripar + c(0, 0, 0, 5))  
  }
  range.h = range(h, na.rm = T)
  range.y = range(y, na.rm = T)
  
  #default valus for ...
  args1=list(xlab="Environment effect",ylab="Variety performance",bty="l",type="n")
  inargs<-list(...)
  args1[names(inargs)]<-inargs
  
  if(!ENVlabel){
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
  h = h[order(h[, chain]), chain]
  
  sorth1 = h[seq(1, length(h), by = 2)]
  sorth2 = h[seq(2, length(h), by = 2)]
  if(ENVlabel){
  axis(side = 1, at = sorth1, labels = names(sorth1), line = 2)
  axis(side = 3, at = sorth2, labels = names(sorth2), line = 1)
  }
  cols = NULL
  pchs = NULL
  for (i in 1:n.VAR) {
    col = i + 1
    pch = 1
    cols = c(cols, col)
    pchs = c(pchs, pch)
    VARi = VARlevels[i]
    wh.i = which(VAR == VARi)
    yhat.i = yhat[wh.i]
    y.i = y[wh.i]
    ENV.i = ENV[wh.i]
    p.i = h[ENV.i]
    order.E.i = order(p.i)
    sortp.i = sort(p.i)
    lines(yhat.i[order.E.i] ~ sortp.i, col = col, type = "l",...)
    points(y.i[order.E.i] ~ sortp.i, col = col, pch = pch,...)
    whmax = which.max(p.i)
    #text(x=min(p.i)+1.05*(max(p.i)-min(p.i)),y=y.i[whmax],labels=IDLi,col=col)
    
  }
  #the intercept denpends on g
  usr <- par("usr")
  clip(range.h[1], range.h[2], range.y[1], range.y[2])
  abline(a = mean(y, na.rm = T), b = 1, lty = 2, col = 1)
  clip(usr[1], usr[2] + (usr[2]-usr[1])* 10, usr[3], usr[4])
  legend(x = usr[2], y = usr[4], legend = c(VARlevels, "slope = 1"), lty = c(rep(1, n.VAR), 2), col = c(cols, 1), bg = "transparent", bty = "n")
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


getIDEL = function(VAR, ENV, VARlevels = NULL, ENVlevels = NULL) {
	VAR = as.character(VAR)
	ENV = as.character(ENV)
	if (is.null(VARlevels)) {
		VARlevels = sort(unique(VAR))
	}
	if (is.null(ENVlevels)) {
		ENVlevels = sort(unique(ENV))
	}
	fVAR = factor(VAR, levels = VARlevels, ordered = T)
	fENV = factor(ENV, levels = ENVlevels, ordered = T)
	IDL = as.numeric(fVAR)
	IDE = as.numeric(fENV)
	out = list()
	out$IDL = IDL
	out$IDE = IDE
	out$VARlevels = VARlevels
	out$ENVlevels = ENVlevels
	out$fVAR = fVAR
	out$fENV = fENV
	return(out)
}



