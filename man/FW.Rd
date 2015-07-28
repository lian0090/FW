\name{FW}
\alias{FW}
\title{
Finlay-Wilkinson Regression
}
\description{
The FW (`Finlay-Wilkinson') function fits Finlay-Wilkinson regression with Ordinary least squares method or Bayesian method (Gibbs Sampler) to continuous traits.
}
\usage{
FW(y, VAR, ENV, method = c("OLS", "Gibbs")[2], saveAt = "", 
    nIter = 5000, burnIn = 3000, thin = 5, df = 5, dfg = 5, dfh = 5, 
    dfb = 5, priorVARe = NULL, priorVARg = NULL, priorVARb = NULL, 
    priorVARh = NULL, A = NULL, nchain = 1, seed = NULL, saveVAR = c(1:2), 
    saveENV = c(1:2)) 
}
\arguments{
  \item{y}{
(numeric, \eqn{n}) The data vector (NAs are alllowed). 
}
  \item{VAR}{
(character, \eqn{n}) is a vector whose elements are identifiers for the varieties which are treated as labels.
}
  \item{ENV}{ 
  (character, \eqn{n}) is a vector whose elements are identifiers for the environments which are treated as labels.
  }
  \item{VARlevels}{(character, \eqn{ng}) is the unique non-duplicated variety identifiers used to specify the order of fitting for the varieties and the order of output for the estimates of \bold{g} and \bold{b}. The order of fitting and output are the same as the order of appearances of variety identifiers in VARlevels; the default values are the unique values of VAR sorted into ascending order by string comparisons. }
  \item{ENVlevels}{(character, \eqn{nh}) is the unique non-duplicated environment identifiers used to specify the order of fitting for the environments and the order of output for the estimates of \bold{h}. The order of fitting and output are the same as the order of appearances of environment identifiers in ENVlevels; the default values are the unique values of ENV sorted into ascending order by string comparison.}
  \item{method}{
  is used to describe what method to use: either "OLS" for ordinary least squares method or "Gibbs" for Gibbs method. The default is "Gibbs".
  }
   \item{A}{ (numeric), is the covariance structure for Normal distirubion of \bold{g} and \bold{b}: N(\bold{g}|\bold{0}, \bold{A}\eqn{\sigma^2_g}),  N(\bold{b}|\bold{0},\bold{A}\eqn{\sigma^2_b}). A may be NULL, this case, \bold{g} and \bold{b} are assumed to be IID}
  \item{saveAt}{ 
  (character) can be used to indicate where to store the samples and to provide a pre-fix to be appended to 'samps.rda' (the name of the file where the samples are stored). By default samples are saved in the current working directory and no pre-fix is added to the file names.   
  }
  \item{nIter, burnIn, thin}{
  (integer) control the number of iterations of the sampler, the number of samples discarded, and the thinning used to compute posterior means, respectively.
  }
  \item{dfe, dfg, dfh, dfb, priorVARe, priorVARg, priorVARb, priorVARh}{
  (numeric) define the hyper parameters of scaled inverse chi-squared distributions for \eqn{\sigma^2_\epsilon}, \eqn{\sigma^2_g}, \eqn{\sigma^2_b} and \eqn{\sigma^2_h}, with dfe, dfg, dfb, dfh degrees of freedom respectively. priorVARe, priorVARg, priorVARb, and  priorVARh defined the prior estimates for  \eqn{\sigma^2_\epsilon}, \eqn{\sigma^2_g}, \eqn{\sigma^2_b} and \eqn{\sigma^2_h}. For a scaled inverse chi-squared variable \eqn{\sigma^2} with degrees of freedom df=\eqn{\nu} and scale parameter \eqn{S^2}, the scale paramter \eqn{S^2} in FW is obtained by satisfying the mode of the density: \eqn{\frac{\nu S^2}{\nu+2}}{\nuS^2/(\nu+2)}. If the prior estimate of \eqn{\sigma^2} is set to be priorVAR=\eqn{\hat{\sigma}^2}, then \eqn{S^2=\hat{\sigma^2}\frac{\nu+2}{\nu}}{S^2=\hat{\sigma}^2(\nu+2)/\nu}. 
  	}
  \item{nchain}{
  (integer) specifies the number of chains for Gibbs Sampler to run. 
  }
  \item{seed}{ 
  (integer vector whose length is equal to the number of chains) is the starting seed for Gibbs Sampler for each chain. If seed=NULL, no seed is set. 
  }
  \item{saveVAR, saveENV}{
  (integer) can be used to specify for which variety or environment the samples of parameter should be saved. For example, saveVAR=c(1,5,10) will save the samples of g and b for variety 1, 5, and 10; saveENV=c(1,5,10) will save the samples of h for environment 1, 5 and 10. By default, only the parameter samples for the first two varieties and first two environments are saved. 
  }
}
\details{
FW implements Gibbs sampler (Gibbs) or Ordinary least square (OLS) method to fit the Finlay-Wilkinson regression model. The basic form of statistical model for both methods is as following:
\deqn{y_{ijk} = \mu  + g_i + h_j + {b_i}{h_j}+\epsilon_{ijk}}{y_ijk = \mu  + g_i + h_j + b_i * h_j + \epsilon_ijk}

\eqn{y_{ijk}}{y_ijk} is the phenotype performance for variety \eqn{i}, environment \eqn{j} and replicate \eqn{k}. \eqn{\mu} is the overall mean, \eqn{g_i} is the effect of variety \eqn{i}, \eqn{h_j} is the effect of environment \eqn{j}, \eqn{b_i}  is the slope effect of variety \eqn{i} on all the environment gradients,\eqn{\epsilon_{ijk}}{\epsilon_ijk} is the normal residual with mean zero and variance \eqn{\sigma^2_\epsilon} . When there are no replicates under variety environment combination, index k can be removed and the above equation reduces to \eqn{y_{ij} = \mu  + {g_i} + {h_j} + {b_i}{h_j}+\epsilon_{ij}}{y_ij = \mu  + g_i + h_j + b_i * h_j +\epsilon_ij}. 

Only the arguments \code{y}, \code{VAR}, \code{ENV} and \code{method} are used for fitting the OLS model. All other arguments are for the Gibbs method.  

Details for model setup can be found in the vignettes in the package. 

}
\value{
A list with estimated posterior means and the parameters used to fit the model.
  \item{y}{response vector used in the call to FW}
  \item{whichNa}{}
  \item{var_g,var_b,var_h,var_e}{the position of the entries in y that were missing}
  \item{VAR}{The identifiers for varieties}
  \item{ENV}{the identifiers for environments}
  \item{VARlevels}{the levels of the varieties (the varieties are fitted in the order of the VARlevels)}
  \item{ENVlevels}{
  levels of the environments (the environments are fitted in the order of the ENVlevels
  }
  \item{mu,g,b,h}{are the parameter estimates (OLS method) or posterior means (Gibbs method) for \eqn{\mu}, \bold{g}, \bold{b}, \bold{h}; mu is a vector; with OLS method, mu has only one element;with Gibbs method, each element of mu represents a different chain. g, b, h are all matrix; with OLS method, g, b, h all have only one column; wiht Gibbs method, each column of g, b, h represents a different chain. 
	}
  \item{yhat}{(matrix) for both methods was the linear predictor obtained from estimated model parameters as \eqn{\hat{y}=\hat{\mu}+\hat{h}+\hat{b}\hat{h}}. For OLS, there is only one column for yhat. For Gibbs, each column of yhat represents a different chain.
  }
  \item{var_e,var_g,var_b,var_h}{are for Gibbs method only, and are the posterior means for \eqn{\sigma^2_\epsilon}, \eqn{\sigma^2_g}, \eqn{\sigma^2_b} and \eqn{\sigma^2_h}.
}
\item{post_yhat}{(matrix) is for Gibbs method only, and is the posterior means of the linear predictor sampled from the corresponding fully conditional density.
}
}
\author{
Lian Lian, Gustavo de los Campos
}
\examples{
\dontrun{
data(wheat.rda)
attach(wheat.Y50)
lm1=FW(y=y,VAR=VAR,ENV=ENV,method="OLS")
lm2=FW(y=y,VAR=VAR,ENV=ENV,method="Gibbs")
}
}
\keyword{models}