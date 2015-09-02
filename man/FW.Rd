\name{FW}
\alias{FW}
\title{
Finlay-Wilkinson Regression
}
\description{
The FW ("Finlay-Wilkinson") function fits the Finlay-Wilkinson regression with Ordinary least squares (OLS) method or Bayesian method (Gibbs Sampler) to continuous traits.
}
\usage{
function (y, VAR, ENV, VARlevels = NULL, ENVlevels = NULL, method = c("OLS", 
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
(character, \eqn{n}) is a vector whose elements are identifiers for the varieties which are treated as labels.
}
  \item{ENV}{ 
  (character, \eqn{n}) is a vector whose elements are identifiers for the environments which are treated as labels.
  }
  \item{VARlevels}{(character, \eqn{qg}) is the unique non-duplicated variety identifiers used to specify the order of fitting for the varieties and the order of output for the estimates of \bold{g} and \bold{b}. The order of fitting and output are the same as the order of appearances of variety identifiers in VARlevels; the default values are the unique values of VAR sorted into ascending order by string comparisons. }
  \item{ENVlevels}{(character, \eqn{qh}) is the unique non-duplicated environment identifiers used to specify the order of fitting for the environments and the order of output for the estimates of \bold{h}. The order of fitting and output are the same as the order of appearances of environment identifiers in ENVlevels; the default values are the unique values of ENV sorted into ascending order by string comparison.}
  \item{method}{
  is used to describe what method to use: either "OLS" for ordinary least squares method or "Gibbs" for Gibbs method. The default is "Gibbs".
  }
   \item{A}{ (numeric), is the covariance structure for Normal distirubion of \bold{g} and \bold{b}: N(\bold{g}|\bold{0}, \bold{A}\eqn{\sigma^2_g}),  N(\bold{b}|\bold{0},\bold{A}\eqn{\sigma^2_b}). A may be NULL, this case, \bold{g} and \bold{b} are assumed to be IID.}
   \item{H}{ (numeric), is the covariance structure for Normal distirubion of \bold{h}: N(\bold{h}|\bold{0}, \bold{H}\eqn{\sigma^2_h}). H may be NULL, this case, \bold{g} and \bold{b} are assumed to be IID.}
  \item{saveAt}{ 
  (character) can be used to indicate where to store the samples and to provide a pre-fix to be appended to 'samps.rda' (the name of the file where the samples are stored). By default samples are saved in the current working directory and no pre-fix is added to the file names.   
  }
  \item{nIter, burnIn, thin}{
  (integer) control the number of iterations of the sampler, the number of samples discarded, and the thinning used to compute posterior means, respectively.
  }
  \item{dfe, dfg, dfh, dfb, priorVARe, priorVARg, priorVARb, priorVARh}{
  (numeric) define the degrees of freedom (df) and prior estimates of variance components (priorVAR) for scaled inverse chi-squared prior distributions:  dfe and priorVARe for \eqn{\sigma^2_\epsilon}, dfg and priorVARg for \eqn{\sigma^2_g}, dfb and priorVARb for \eqn{\sigma^2_b}, dfh and priorVARh for \eqn{\sigma^2_h}. For a scaled inverse chi-squared variable \eqn{\sigma^2} with degrees of freedom df and scale parameter \eqn{S^2}, the scale paramter \eqn{S^2} in FW is obtained by satisfying the mode of the density: \eqn{\frac{\nu S^2}{df+2}}{df*S^2/(df+2)}. If the prior estimate of \eqn{\sigma^2} is set to be priorVAR, then \eqn{S^2=priorVAR\frac{\nu+2}{\nu}}{S^2=priorVAR(df+2)/df}. 
  	}
  \item{nchain}{
  (integer) specifies the number of chains for Gibbs Sampler to run. 
  }
  \item{seed}{ 
  (integer vector whose length is equal to the number of chains) is the starting seed for Gibbs Sampler for each chain. If seed=NULL, no seed is set. 
  }
  \item{inits}{ (list of named lists) specifies the initial values for Gibbs sampler. The parameters taken by R (the corresponding math representation in in parentheses) are: mu (\eqn{\mu}), g (\bold{g}), b(\bold{b}), h(\bold{h}), var_e(\eqn{\sigma^2_\epsilon}), var_g(\eqn{\sigma^2_g}), var_b(\eqn{\sigma^2_b}) and var_h(\eqn{\sigma^2_h}). These values must be supplied in a list as inits=list(inits1=list(mu=,g=,b=,h=,var_e=,var_g=,var_b=,var_h=),inits2=list(...)). The default initial values for chain 1 are setting mu, g, b, h all to 0, and setting var_e as 0.5 var_y ,var_g as 0.25 var_y, var_b and var_h as 0.5 sd_y, where var_y and sd_y are the total phenotypic variance and standard deviation. The default initial values for all other chains are set by random sampling: mu from N(0, 0.5 var_y), g from N(0, 0.25 var_y ), b and h from N(0, 0.5 sd_y ), var_e, var_g , var_b and var_h  all from uniform distributions with the boundaries set as half and twice the values used in chain 1: var_e from  U(0.5, 2)0.5var_y , var_g from  U(0.5, 2)0.25var_y, var_b and var_h from U(0.5, 2)0.5sd_y.
  }
  \item{saveVAR, saveENV}{
  (integer) can be used to specify for which variety or environment the samples of parameter should be saved. For example, saveVAR=c(1,5,10) will save the samples of g and b for variety 1, 5, and 10; saveENV=c(1,5,10) will save the samples of h for environment 1, 5 and 10. By default, only the parameter samples for the first two varieties and first two environments are saved. 
  }
}
\details{
\bold{Model Specification}.
FW implements Gibbs sampler (Gibbs) or Ordinary least square (OLS) method to fit the Finlay-Wilkinson regression model. The basic form of statistical model for both OLS and Gibbs methods is the same: 
\deqn{y_{ijk} = \mu  + g_i + h_j + {b_i}{h_j}+\epsilon_{ijk} [1]}{y_ijk = \mu  + g_i + h_j + b_i * h_j + \epsilon_ijk [1]} 

\eqn{y_{ijk}}{y_ijk} is the phenotype performance for variety \eqn{i}, environment \eqn{j} and replicate \eqn{k}. \eqn{\mu} is the overall mean, \eqn{g_i} is the effect of variety \eqn{i}, \eqn{h_j} is the effect of environment \eqn{j}, \eqn{b_i}  is the slope effect of variety \eqn{i} on all the environment gradients,\eqn{\epsilon_{ijk}}{\epsilon_ijk} is the normal residual with mean zero and variance \eqn{\sigma^2_\epsilon} . When there are no replicates under variety environment combination, index k can be removed and the above equation reduces to \eqn{y_{ij} = \mu  + {g_i} + {h_j} + {b_i}{h_j}+\epsilon_{ij}}{y_ij = \mu  + g_i + h_j + b_i * h_j +\epsilon_ij}. 

\bold{Posterior distribution,likelihood and prior distributions}. Inferences in a Bayesian model are based on the posterior distribution of unknowns given the data: p(\bold{θ}│\bold{y})~p(\bold{y}|\bold{θ})p(\bold{θ}), where \bold{θ} represents the collection of the unknowns: μ, \bold{g}, \bold{b}, \bold{h}, and \eqn{\sigma^2_\epsilon}, p(\bold{y}|\bold{θ}) is the likelihood. Accoridng to Eqation [1] and assuming IID normal residuals, we have (\bold{y}|\bold{θ})=\eqn{\prod_ijk} N(\eqn{\mu}+g_i + h_j + b_i h_j,\eqn{\sigma^2_\epsilon}).

In the FW package, the prior density is assumed to have the following form: p(\bold{θ})=p(\eqn{\mu})p(\eqn{\sigma^2_\epsilon})p(\bold{b})p(\bold{h})p(\bold{b}). The intercept is assigned a flat prior and the residual variance \eqn{\sigma^2_\epsilon} is assigned a scaled-inverse chi-squared distribution: \eqn{\sigma^2_\epsilon} ~ (\eqn{\nu_\epsilon},\eqn{S^2_\epsilon}), with degrees of freedom \eqn{\nu_\epsilon} (>0) and scale parameter \eqn{S^2_\epsilon} (>0). The prior for \bold{h}, \bold{g} and \bold{b} are all multivariate Normal: \bold{h}~N(\bold{0},\bold{H}\eqn{\sigma^2_h}), \bold{g}~N(\bold{0},\bold{A}\eqn{\sigma^2_g}), \bold{b}~N(\bold{0},\bold{A}\eqn{\sigma^2_b}). Since σ_h^2, σ_g^2 and σ_b^2 are unknown, they are also assigned scaled-inverse chi-squared distributions whose shape are controlled by a set of hyper parameters.




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
  \item{mu,g,b,h}{are the parameter estimates (OLS method) or estimates of posterior means (Gibbs method) for \eqn{\mu}, \bold{g}, \bold{b}, \bold{h}; mu is a vector; with OLS method, mu has only one element;with Gibbs method, each element of mu represents a different chain. g, b, h are all matrix; with OLS method, g, b, h all have only one column; wiht Gibbs method, each column of g, b, h provides estimates derived from one MCMC chain. 
	}
  \item{yhat}{(matrix) is the parameter estimates (OLS method) or estimates of posterior means (Gibbs method) of the predictor \eqn{\hat{y}}: \eqn{\hat{y_ijk}=\hat{\mu}+\hat{g_i}+\hat{h_j}+\hat{b_i}\hat{h_j}}. For OLS, there is only one column for yhat. For Gibbs, each column of yhat represents a different chain.
  }
  \item{var_e,var_g,var_b,var_h}{are for Gibbs method only, and are the posterior means for \eqn{\sigma^2_\epsilon}, \eqn{\sigma^2_g}, \eqn{\sigma^2_b} and \eqn{\sigma^2_h}.
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