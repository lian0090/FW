\name{wheat.G}
\alias{wheat.G}
\alias{wheat.G50}
\title{Pedigree info for the wheat dataset}
\description{
wheat.G (matrix, 599 x 599),  was calculated from DArT markers: marker genotyped were pruned based LD with   threshold of 0.2, and implemented by plink (Purcell et al., 2007) with the options "--indep-pairwise 100 10 0.2". The pruned markers were then centered and standardized to a matrix \eqn{\tilde{X}}, G matrix was calculated as: \eqn{\tilde{X}\tilde{X}'/p} , where p is the number of columns in \eqn{\tilde{X}}. 

wheat.G50 is a subset of wheat.G, and is the genomic matrix of the 50 randomly selected wheat lines in \code{\link{wheat.Y50}}. 

}

\source{
  International Maize and Wheat Improvement Center (CIMMYT), Mexico.
}

\references{
McLaren, C. G., R. Bruskiewich, A.M. Portugal, and A.B. Cosico. 2005. The International Rice Information System. 
A platform for meta-analysis of rice crop data. \emph{Plant Physiology} \bold{139}: 637-642.
}
\keyword{datasets}
