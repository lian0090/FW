# FW
Bayesian method and ordinary least square method for Finlay-Wilkinson Regression. 

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

# Note
The FW regression was fitting the model y=mu+g+(b+1)h+e
In Ordianary Least Square method, a linear regression model is fitted within each line/variety to estimate the genetic effect g and the slope (b+1) where mu is set to 0 for this within line/variety linear regression and g is the intercept of each within line regression. However, in the Bayesian method, the whole data set was used to fit the model, and the g is estimated as a random effect with zero mean. Therefore, the estimated g from ordinary linear regression is generally much larger than the g estimated from Bayesian method due to fact that the estimated g from ordinary linear regression contains an overall mean.  


Detailed implementation can be found in the vignettes of the package. 


