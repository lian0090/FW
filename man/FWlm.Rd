\name{FWlm}
\alias{FWlm}
\title{
Ordinary Least Squares to Perform Finley-Wilkinson Regression: environment means are the environment gradients
}
\description{
Ordinary Least Squares to Perform Finley-Wilkinson Regression: environment means are the environment gradients
}
\usage{
FWlm(y, IDL, IDE)
}
\arguments{
  \item{y}{
Response variable (observations).
}
  \item{IDL}{
Integer values that relate the genotypes (lines) to observations.
}
  \item{IDE}{
Integer values that relate the Environment to the observations.
}
}
\details{
y=g[IDL]+h[IDE]+b[IDL]*h[IDE]

The environment gradient h is the environment means, g and b are obtained by regression of y on h. 

}
\value{
a list with estimated values for g, b, h. 
}
\references{
Finlay, K. W. and G. N. Wilkinson. 1963. The analysis of adaptation in a plant-breeding programme. \emph{Crop and Pasture Science}, \bold{14}: 742-754.

}
\author{
Lian Lian
}


\examples{
data(FWdat)
dat=FWdat$dat
FWlm(y=dat$y,dat$IDL,dat$IDE)
}
\keyword{models}

