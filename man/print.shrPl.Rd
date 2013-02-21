\name{print.shrPl}
\alias{print.shrPl}

\title{
Print method for \code{shrPl} objects 
}

\usage{
\method{print}{shrPl}(x, conf.int = 0.95, digits = 4,
 pvalDigits = 4, eps = 1e-04, ...)
}
\description{Print a summary of a fitted illness-death model using the penalized likelihood approach.}
\arguments{
  \item{x}{
a \code{shrPl} object, i.e., the result of a call to the \code{\link{shr}} function with \code{hazard}="Splines".
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
 \code{\link{summary.shrPl}},  
 \code{\link{plot.shrPl}}
}
  
\author{
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>

}

\examples{
\dontrun{
# a penalized survival model 
data(test)
fit.su <- shr(Hist(time=list(l,r),id)~cov,data=test,hazard="Splines") 
print(fit.su) 
}
}

\keyword{methods}
