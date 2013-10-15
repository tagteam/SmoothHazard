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
a \code{idmWeib} object, i.e., the result of a call to the \code{\link{idm}} function with \code{intensities}="Weib".
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
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>

}

\examples{
\dontrun{
# Weibull illness-death model 
data(Paq1000)

fit <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
formula01=Hist(time=list(l,r),event=dementia)~certif,data=Paq1000) 

}
}

\keyword{methods}
