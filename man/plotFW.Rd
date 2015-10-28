\name{plot.FW}
\alias{plot.FW}
\title{
plot the FW object 
}
\description{
the plot.FW function plot the fitted object from calling the function \code{\link{FW}}.
}

\usage{
plot.FW(FWobj,plotVAR=NULL,chain=1,...)
}
\arguments{
\item{FWojb}{(class 'FW') returned from calling \code{\link{FW}} function.}
\item{plotVAR}{(character vector or integer vector) If plotVAR is a character vector, it specifies the variety identifiers for which the fitted values should be plotted. If plotVar is an integer vector, it specifies  the indexes of the varieties for which the fitted values should be plotted (for example, plotVAR=c(1:5) will plot the first five varieties). The default values is NULL and all the varieties will be plotted.
}
\item{chain}{(integer), only used for Gibbs method, when there are fitted values from multiple chains. This specifies from which chain the fitted value is to be plotted
}
\item{...}{Arguments to be passed for plot function (e.g., xlab, ylab, main, ...).}
}

\details{
The x-axis is the estimated values of environment effects. The y-axis is the variety performance. The colored lines are the fitted values against environment gradients. The colored circles are the observed cell means for the variety environment combinations. 
}
\author{
Lian Lian, Gustavo de los Campos
}
\examples{
\dontrun{
data(wheat)
attach(wheat.Y)
lm2=FW(y=y[1:100],VAR=VAR[1:100],ENV=ENV[1:100])
plot(lm2)
plot(lm2,plotVAR=c("1081265","1101307", "1295736", "13302", "1343502"))
}
}