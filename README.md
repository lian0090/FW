# FW
Gibbs Sampler for Finlay-Wilkinson Regression

# Install
```R
library(devtools)
install_github("lian0090/FW")
```
Note, install_github by default does not build vignettes. To be able to directly view vignettes from R, you need to do the folowing installation isntead. 
`install_github("lian0090/FW",build_vignettes=T)`

# Basic usage
```R
library(FW)
data(wheat)
attach(wheat.Y)
lm1=FW(y=y,VAR=VAR,ENV=ENV)
plot(lm1)
```


