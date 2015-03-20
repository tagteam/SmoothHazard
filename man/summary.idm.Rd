\name{summary.idm}
\alias{summary.idm}
\title{Summary of a fitted illness-death model}
\usage{
\method{summary}{idm}(object, conf.int = 0.95, digits = 4, pvalDigits = 4,
  eps = 0.0001, ...)
}
\arguments{
  \item{object}{a \code{idmSplines} object, i.e., the
  result of a call to the \code{\link{idm}} function with
  \code{intensities}="Splines".}

  \item{conf.int}{confiance level.}

  \item{digits}{number of digits to print.}

  \item{pvalDigits}{number of digits to print for
  p-values.}

  \item{eps}{convergence criterion used for p-values.}

  \item{\dots}{other unusued arguments.}
}
\description{
Summarize the event history data of an illness-death
regression model and show regression coefficients for
transition intensities
}
\examples{
\dontrun{
library(prodlim)
data(Paq1000)
fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=e)~certif,
		formula01=Hist(time=list(l,r),event=dementia)~certif,
                formula12=~1,
                method="Splines",
		data=Paq1000)
summary(fit.splines)
}
}
\author{
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>
}
\seealso{
\code{\link{idm}}, \code{\link{print.idm}},
\code{\link{plot.idm}}
}
\keyword{methods}

