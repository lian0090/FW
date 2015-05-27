\name{FWdat}
\alias{FWdat}
\docType{data}
\title{
A simulated dataset for Finley-Wilkinson regression.
}
\description{
A simulated dataset for Finley-Wilkinson regression.
}
\usage{data(FWdat)}
\format{
  The format is:
List of 2
 $ dat          a data.frame with 3 variables:
  ..$ y  : observations
  ..$ IDL: integer vector that relates gentoypes (Lines) to observations
  ..$ IDE: integer vector that relates Environments to observations
 $ realizedValue: the realized values of the simulated parameters.
  ..$ mu   : overall mean
  ..$ g    : genotype main effects
  ..$ h    : environment main effects
  ..$ b    : slope for genotypes in Finley-Wilkinson regression
  ..$ var_g: g ~ N(0,var_g)
  ..$ var_h: h ~ N(0,var_h)
  ..$ var_b: b ~ N(0,var_b) 
  ..$ var_e: e ~ N(0,var_e)
}
\details{
y=mu+g[IDL]+b[IDL]*h[IDE]+h[IDE]+e
simulate g from N(0, var_g), b from N(0,var_b), h from N(var_h) and e from N(0,var_e). var_g=var(g), var_h=var(h), var_b=var(b) and var_e=var(e) 
}
\keyword{datasets}
