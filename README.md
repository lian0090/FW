# FW
Bayesian method and ordinary least square methods for Finlay-Wilkinson Regression. 

# Install
```R
library(devtools)
install_github("lian0090/FW")
```
Note, install_github by default does not build vignettes. To be able to directly view vignettes from R, you need to do the following installation instead. 
`install_github("lian0090/FW",build_vignettes=T)`

# Basic usage
```R
library(FW)
data(wheat)
attach(wheat.Y)
lm1=FW(y=y,VAR=VAR,ENV=ENV)
plot(lm1)
```

Detailed implementation can be found in the vignettes of the package. 


