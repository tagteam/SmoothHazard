\name{summary.shr}
\alias{summary.shr}
\title{Summary of a fitted survival model using a penalized likelihood approach}
\usage{
\method{summary}{shr}(object, conf.int = 0.95, digits = 4, pvalDigits = 4,
  eps = 0.0001, ...)
}
\arguments{
  \item{object}{a \code{shr} object, i.e., the result of a
  call to the \code{\link{shr}} function.}

  \item{conf.int}{The confidence level.}

  \item{digits}{number of digits to print.}

  \item{pvalDigits}{number of digits to print for
  p-values.}

  \item{eps}{convergence criterion used for p-values.}

  \item{\dots}{other unusued arguments.}
}
\description{
Print a short summary of a fitted illness-death model using
the penalized likelihood approach.
}
\examples{
\dontrun{
# a penalized survival model
data(testdata)
fit.su <- shr(Hist(time=list(l,r),id)~cov,data=testdata,method="Splines")
summary(fit.su)

# Weibull survival model
data(testdata)
fit.su <- shr(Hist(time=list(l,r),id)~cov,data=testdata)
summary(fit.su)
}
}
\author{
Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
}
\seealso{
\code{\link{shr}}, \code{\link{print.shr}},
\code{\link{plot.shr}}
}
\keyword{methods}

