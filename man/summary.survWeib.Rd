\name{summary.survWeib}
\alias{summary.survWeib}

\title{
Summary of a fitted survival model using a Weibull approach
}

\usage{
\method{summary}{survWeib}(object, conf.int = 0.95, digits = 4,
 pvalDigits = 4, eps = 1e-04, ...)
}
\description{Print a short summary of a fitted illness-death model using the Weibull approach.}
\arguments{
  \item{object}{
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
 \code{\link{surv}}, 
 \code{\link{print.survWeib}},  
 \code{\link{plot.survWeib}}
}
  
\author{
J. Pierre
}

\examples{

# Weibull survival model 
data(test)
fit.su <- surv(Hist(time=list(l,r),id)~cov,data=test) 
summary(fit.su) 
}


\keyword{methods}
