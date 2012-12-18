\name{print.survWeib}
\alias{print.survWeib}

\title{
Print method for \code{survWeib} objects 
}

\usage{
\method{print}{survWeib}(x, conf.int = 0.95, digits = 4,
 pvalDigits = 4, eps = 1e-04, ...)
}
\description{Print a summary of a fitted illness-death model using the penalized likelihood approach.}
\arguments{
  \item{x}{
a \code{survWeib} object, i.e., the result of a call to the \code{\link{surv}} function with \code{hazard}="Weib".
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
 \code{\link{summary.survWeib}},  
 \code{\link{plot.survWeib}}
}
  
\author{
J. Pierre
}

\examples{
\dontrun{
#  Weibull survival model 
data(test)
fit.su <- surv(Hist(time=list(l,r),id)~cov,data=test) 
print(fit.su) 
}
}

\keyword{methods}
