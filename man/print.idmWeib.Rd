\name{print.idmWeib}
\alias{print.idmWeib}

\title{
Print method for \code{idmWeib} objects 
}

\usage{
\method{print}{idmWeib}(x, conf.int = 0.95, digits = 4,
 pvalDigits = 4, eps = 1e-04, ...)
}
\description{Print a summary of a fitted illness-death model using the Weibull approach.}
\arguments{
  \item{x}{
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
 \code{\link{summary.idmWeib}},  
 \code{\link{plot.idmWeib}}
}
  
\author{
J. Pierre
}

\examples{
\dontrun{
# Weibull survival model ***/
data(test)
fit.su <- shr(Hist(time=list(l,r),id)~cov,data=test) 
print(fit.su) 
}
}

\keyword{methods}
