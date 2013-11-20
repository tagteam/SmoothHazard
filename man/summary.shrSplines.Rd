\name{summary.shrSplines}
\alias{summary.shrSplines}

\title{
Summary of a fitted survival model using a penalized likelihood approach
}

\usage{
\method{summary}{shrSplines}(object, conf.int = 0.95, digits = 4,
 pvalDigits = 4, eps = 1e-04, ...)
}
\description{Print a short summary of a fitted illness-death model using the penalized likelihood approach.}
\arguments{
  \item{object}{
a \code{shrSplines} object, i.e., the result of a call to the \code{\link{shr}} function with \code{hazard}="Splines".
}
  \item{conf.int}{
confiance level.
}
  \item{digits}{
number of digits to print.
}
  \item{pvalDigits}{
number of digits to print for p-values.
}
  \item{eps}{
convergence criterion used for p-values.
}
  \item{\dots}{
other unusued arguments.
}
}

\seealso{
 \code{\link{shr}}, 
 \code{\link{print.shrSplines}},  
 \code{\link{plot.shrSplines}}
}
  
\author{
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>

}

\examples{
\dontrun{
# a penalized survival model 
data(testdata)
fit.su <- shr(Hist(time=list(l,r),id)~cov,data=testdata,hazard="Splines") 
summary(fit.su) 
}
}

\keyword{methods}
