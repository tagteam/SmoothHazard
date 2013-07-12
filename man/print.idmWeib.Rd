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
    Class \code{idmWeib} object, i.e. the result of a call to the \code{\link{idm}} function with \code{intensities}="Weib".
  }
  \item{conf.int}{
    Confidence level.
  }
  \item{digits}{
    Number of digits to print.
  }
  \item{pvalDigits}{
    Number of digits to print for p-values.
  }
  \item{eps}{
    Passed to \code{format.pval}.
  }
  \item{\dots}{
    Not used.
  }
}

\seealso{
 \code{\link{summary.idmWeib}},  
 \code{\link{plot.idmWeib}}
}
  
\author{
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>

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
