\name{FW}
\alias{FW}
\title{
An R package for Genomic/Pedigree and Spatial analysis using the Finlay-Wilkinson Regression
}
\description{
The FW ("Finlay-Wilkinson") function fits the Finlay-Wilkinson regression with Ordinary least squares method or Bayesian method to continuous traits. Genomic informationa and Spatial correlations can be included in the model. 
}
\usage{
FW(y, VAR, ENV, VARlevels = NULL, ENVlevels = NULL, method = c("OLS", 
    "Gibbs")[2], A = NULL, H = NULL, saveAt = "", nIter = 5000, 
    burnIn = 3000, thin = 5, dfe = 5, dfg = 5, dfh = 5, dfb = 5, 
    priorVARe = NULL, priorVARg = NULL, priorVARb = NULL, priorVARh = NULL, 
    nchain = 1, seed = NULL, inits = NULL, saveVAR = c(1:2), 
    saveENV = c(1:2)) 
}
\arguments{
Only the arguments \code{y}, \code{VAR}, \code{ENV}, \code{VARlevels}, \code{ENVlevels} and \code{method} are used for fitting the OLS model. All other arguments are for the Gibbs method.  
  \item{y}{
(numeric, \eqn{n}) The data vector (NAs are alllowed). 
}
  \item{VAR}{
(character, \eqn{n}) is a vector whose elements are identifiers for the varieties which are treated as labels (NAs are not allowed).
}
  \item{ENV}{ 
  (character, \eqn{n}) is a vector whose elements are identifiers for the environments which are treated as labels (NAs are not allowed).
  }
  \item{VARlevels}{(character, \eqn{qg}) unique levels of variety identifiers used to specify the order of fitting for the varieties and the order of output for the estimates of g and b. The default values are the unique values of VAR sorted into ascending order by string comparisons. }
  \item{ENVlevels}{(character, \eqn{qh}) unique levels of environment identifiers used to specify the order of fitting for the environments and the order of output for the estimates of h. The default values are the unique values of ENV sorted into ascending order by string comparison.}
  \item{method}{
  Describes what method to use: either "OLS" for ordinary least squares method or "Gibbs" for Bayesian method using Gibbs Sampler. The default is "Gibbs".
  }
   \item{A}{ (numeric, \eqn{qg*qg}), is the covariance structure for Normal distirubion of g and b: N(0, A\eqn{\sigma^2_g}),  N(b|0,A\eqn{\sigma^2_b}). In case A is NULL, A is automatically set to be Identiy matrix}
   \item{H}{ (numeric,\eqn{qh*qh}), is the covariance structure for Normal distirubion of h: N(0, H\eqn{\sigma^2_h}). In case H is NULL, H will be automatically set to be Identiy matrix }
  \item{saveAt}{ 
  (character) can be used to indicate where to store the samples and to provide a pre-fix to be appended to 'samps.rda' (the name of the file where the samples are stored). By default samples are saved in the current working directory and no pre-fix is added to the file names.   
  }
  \item{nIter, burnIn, thin}{
  (integer) control the number of iterations of the sampler, the number of samples discarded, and the thinning used to compute posterior means.
  }
  \item{dfe, dfg, dfh, dfb, priorVARe, priorVARg, priorVARb, priorVARh}{
  (numeric) define the degrees of freedom (df) and prior estimates of variance components (priorVAR) for scaled inverse chi-squared distributions. The scale paramter \eqn{S^2} is set to be priorVAR(df+2)/df. Details can be found in the vignettes of the software. 
  	}
  \item{nchain}{
  (integer) specifies the number of chains for Gibbs Sampler to run. 
  }
  \item{seed}{ 
  (integer vector whose length is equal to the number of chains) is the starting seed for Gibbs Sampler for each chain. If seed=NULL, no seed is set. 
  }
  \item{inits}{ (list of named lists) specifies the initial values for Gibbs sampler. Example usage is inits=list(inits1=list(mu=,g=,b=,h=,var_e=,var_g=,var_b=,var_h=),inits2=list(...)). If NULL, default values will be set. Details can be seen in the vignettes of the software. 
  }
  \item{saveVAR, saveENV}{
  (integer) can be used to specify for which variety or environment the samples of parameter should be saved. For example, saveVAR=c(1,5,10) will save the samples of g and b for variety 1, 5, and 10. By default, only the parameter samples for the first two varieties and first two environments are saved. 
  }
}
\details{
\bold{Model Specification}.
FW implements Gibbs sampler (Gibbs) or Ordinary least square (OLS) method to fit the Finlay-Wilkinson regression model. The basic form of statistical model for both OLS and Gibbs methods is the same: 
\deqn{y = \mu  + g + h + bh+\epsilon [1]}{y= \mu  + g + h + b * h + \epsilon [1]} 

\eqn{y}{y} is the phenotype performance, \eqn{\mu} is the overall mean, \eqn{g} is the main variety effect, \eqn{h} is the main environment effect, \eqn{b} is the variety slope on the environment gradients, \eqn{\epsilon}{\epsilon} is the normal residual with mean zero and variance \eqn{\sigma^2_\epsilon}.


Further details for Model setup and implementation can be found in the vignettes. 



}
\value{
A list with estimated posterior means and the parameters used to fit the model.
  \item{y}{response vector used in the call to FW}
  \item{whichNa}{the position of the entries in y that were missing}
  \item{VAR}{The identifiers for varieties}
  \item{ENV}{the identifiers for environments}
  \item{VARlevels}{the levels of the varieties (the varieties are fitted in the order of the VARlevels)}
  \item{ENVlevels}{
  levels of the environments (the environments are fitted in the order of the ENVlevels
  }
  \item{mu,g,b,h}{are the estimated parameters. mu is a vector: with OLS method, mu has only one element; with Gibbs method, each element of mu is an estimate from a different chain. g, b, h are all matrix: with OLS method, g, b, h all have only one column; with Gibbs method, each column of g, b, h provides estimates derived from one MCMC chain. 
	}
  \item{yhat}{(matrix) is the estimates of the predictor: \eqn{yhat = mu+ g + h + b*h}. For OLS, there is only one column for yhat. For Gibbs method, each column of yhat corresponds to a different chain.
  }
  \item{SD.mu,SD.g,SD.b,SD.h,SD.yhat}{only for Gibbs method, the posterior standard deviation for mu, g, b, h and yhat}
  \item{var_e,var_g,var_b,var_h}{only for Gibbs method, the posterior means for \eqn{\sigma^2_\epsilon}, \eqn{\sigma^2_g}, \eqn{\sigma^2_b}, \eqn{\sigma^2_h}
  }
  \item{SD.var_e,SD.var_g,SD.var_b,SD.var_h}{only for Gibbs method, the posterior standard deviation for \eqn{\sigma^2_\epsilon}, \eqn{\sigma^2_g}, \eqn{\sigma^2_b}, \eqn{\sigma^2_h}
}
}
\author{
Lian Lian, Gustavo de los Campos
}
\examples{
data(wheat)
attach(wheat.Y)
lm1=FW(y=y[1:100],VAR=VAR[1:100],ENV=ENV[1:100],method="OLS")
lm2=FW(y=y[1:100],VAR=VAR[1:100],ENV=ENV[1:100],method="Gibbs")
}
\keyword{models}