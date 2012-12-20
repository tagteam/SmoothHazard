\name{summary.idmWeib}
\alias{summary.idmWeib}

\title{
Summary of a fitted illness-death model using a Weibull approach
}

\usage{
\method{summary}{idmWeib}(object, conf.int = 0.95, digits = 4,
 pvalDigits = 4, eps = 1e-04, ...)
}
\description{Print a short summary of a fitted Weibull model using the penalized likelihood approach.}
\arguments{
  \item{object}{
a \code{idmWeib} object, i.e., the result of a call to the \code{\link{idm}} function with \code{hazard}="Weib".
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
 \code{\link{idm}}, 
 \code{\link{print.idmWeib}},  
 \code{\link{plot.idmWeib}}
}
  
\author{
J. Pierre
}

\examples{
\dontrun{
# Weibull survival model 
data(test)
fit.su <- shr(Hist(time=list(l,r),id)~cov,data=test) 
fit.su 
}
}

\keyword{methods}
