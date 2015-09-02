cd ~/Dropbox/work/arvalisGxE/manuscript/FW/src
R CMD SHLIB sample_beta.c C_GibbsFW.c -o C_GibbsFW.so
FWdir="~/Dropbox/work/arvalisGxE/manuscript/FW"
setwd(file.path(FWdir,"src"))
dyn.load("C_GibbsFW.so")
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
       if(trace) cat(nm,":")
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
 }

sourceDir(file.path(FWdir,"R"))
library(coda)

load("~/Dropbox/work/arvalisGxE/wheat/wheatAll.rda")

attach(wheat.Y50)
lm1=FW(y,VAR,ENV)

H=diag(1,4)
colnames(H)=rownames(H)=unique(ENV)
lm1=FW(y,VAR,ENV,A=wheat.G)
lm2=FW(y,VAR,ENV,A=wheat.G,H=H)
cor(lm1$y,lm1$yhat)
cor(lm2$y,lm2$yhat)

whichNa=which(is.na(yNA))
lm1=FW(yNA,VAR,ENV,A=wheat.G,seed=2)
cor(lm1$yhat[whichNa],y[whichNa])

lm2=FW(yNA,VAR,ENV,A=wheat.G,H=H,seed=2)
cor(lm2$yhat[whichNa],y[whichNa])

lm3=FW(yNA,VAR,ENV,method="OLS")
cor(lm3$yhat[whichNa],y[whichNa])

load("samps.rda")
plot(samps,ask=T,density=F)