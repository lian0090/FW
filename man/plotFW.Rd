\name{plot.FW}
\alias{plot.FW}
\title{
plot the FW object 
}
\description{
the plot.FW function plot the fitted object from calling the function \code{\link{FW}}.
}

\usage{
plot.FW(FWobj,plotVAR=NULL,main=NULL,chain=1)
}
\arguments{
\item{FWojb}{(class 'FW') returned from calling \code{\link{FW}} function.}
\item{plotVAR}{(character vector) the variety identifiers for which the fitted values should be plotted. The plotVAR can be NULL, and all the varieties will be plotted.
}
\item{main}{an overall title for the plot}
\item{chain}{(integer), only used for Gibbs method, when there are fitted values from multiple chains. This specifies from which chain the fitted value is to be plotted}
}
\details{
The x-axis is the estimated values of environment effects. The colored lines are the fitted values against environment gradients. The colored circles are the observed cell means for the variety environment combinations. 
}
\author{
Lian Lian, Gustavo de los Campos
}
\examples{
\dontrun{
data(wheat)
attach(wheat)
lm2=FW(y=y,VAR=VAR,ENV=ENV,method="Gibbs",nIter=200, burnIn=100)
plot(lm2)
plot(lm2,plotVAR=c("1081265","1101307", "1295736", "13302", "1343502"))
}
}