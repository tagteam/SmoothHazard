\name{print.idmSplines}
\alias{print.idmSplines}

\title{
Print method for \code{idmSplines} objects 
}

\usage{
\method{print}{idmSplines}(x, conf.int = 0.95, digits = 4, 
pvalDigits = 4, eps = 1e-04, ...)
}
\description{Print a summary of a fitted illness-death model using the penalized likelihood approach.}
\arguments{
  \item{x}{
a \code{idmSplines} object, i.e., the result of a call to the \code{\link{idm}} function with \code{intensities}="Splines".
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
 \code{\link{summary.idmSplines}},  
 \code{\link{plot.idmSplines}}
}
  
\author{
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>

}

\examples{
\dontrun{
data(Paq1000)
fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=t0)~certif,
		formula01=Hist(time=list(l,r),event=dementia)~certif,
                formula12=~1,
                intensities="Splines",
		data=Paq1000)
print(fit.splines)
}
}

\keyword{methods}
