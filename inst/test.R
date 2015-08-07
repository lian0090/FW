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

lm1=FW(y,VAR,ENV,A=wheat.G)
cor(lm1$y,lm1$yhat)

load("samps.rda")
plot(samps,ask=T,density=F)