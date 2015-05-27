\name{GibbsFW}
\alias{GibbsFW}
\title{
Gibbs Sampler to perform Finley-Wilkinson regression
}
\description{
Gibbs Sampler to perform Finley-Wilkinson regression
}
\usage{
GibbsFW(y = y, IDL, IDE ,savedir=".", nIter = 1000, burnIn = 500, thin = 1, df = 5, dfg = 5, dfh = 5, dfb = 5, S = NULL, Sg = NULL, Sb = NULL, Sh = NULL, A = NULL,inits=NULL,seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
(numeric, \eqn{n}{n}) The data vector (NAs not alllowed). 
}
  \item{IDL}{
(integer, \eqn{n}{n}) The vector linking observations to genotypes (lines). IDL must be an integer taking values from 1 to \eqn{ng}{ng}; IDL[i]=\eqn{q}{q} indicates that the ith observation in \eqn{\boldsymbol y}{y} belongs to genotype (line) \eqn{q}{q}. }
  \item{IDE}{ (integer, \eqn{n}{n}) The vector linking observations to Environment.  IDE must be an integer taking values from 1 to \eqn{nh}{nh}; IDE[i]=\eqn{q}{q} indicates that the ith observation in \eqn{\boldsymbol y}{y} belongs to Environment \eqn{q}{q}.}
  \item{savedir}{directory to save the samples}
  \item{nIter, burnIn, thin}{(integer) the number of iterations, burn-in and thinning.}
  \item{df,dfg,dfb,dfh,S,Sg,Sb,Sh}{The scale parameter for the scaled inverse-chi squared prior assigned to the variance (df,S for var_e), (dfg,Sg for var_g), (dfb,Sb for var_b) and (dfh,Sh for var_h)
  	}
  \item{inits}{Initial values for mu, g, b, h, var_e, var_g,var_b, var_h.Must be supplied in a list as inits=list(inits1=list(mu=,g=,....),inits2=).The default is to set two sets of initial values}
  \item{A}{ a (co)variance structure (numeric, \eqn{ng \times ng}{ng x ng}) between effects of the genotypes (lines).  whose (co)variance function is in the qth row (column) of \eqn{\boldsymbol A}{A}.
            A may be NULL, this case, lines are assumed to be independent.}
  \item{seed}{Random seed to generate inital values for default chain 2 and to do Gibbs Sampling}          
}

\value{
	return a list of posteria means for the parameters: mu, var_g, var_b, var_h, var_e, g, b, and h.
	The samples are saved under the working directory as samps.rda. samps is a mcmc.list type that is ready to be used by the coda package. 
}
\author{
Lian Lian
}
\examples{
	data(FWdat)
	dat=FWdat$dat
	ans=GibbsFW(y=dat$y,IDL=as.integer(as.character(dat$IDL)),IDE=as.integer(as.character(dat$IDE)),savedir=".",nIter=1000,burnIn=500)
}
\keyword{models}