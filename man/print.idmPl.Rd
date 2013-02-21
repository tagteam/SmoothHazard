\name{print.idmPl}
\alias{print.idmPl}

\title{
Print method for \code{idmPl} objects 
}

\usage{
\method{print}{idmPl}(x, conf.int = 0.95, digits = 4, 
pvalDigits = 4, eps = 1e-04, ...)
}
\description{Print a summary of a fitted illness-death model using the penalized likelihood approach.}
\arguments{
  \item{x}{
a \code{idmPl} object, i.e., the result of a call to the \code{\link{idm}} function with \code{hazard}="Splines".
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
 \code{\link{summary.idmPl}},  
 \code{\link{plot.idmPl}}
}
  
\author{
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>

}

\examples{
\dontrun{
data(Paq1000)
d <- Paq1000
names(d) <- c("dementia","mort","entry","L","R","time","certif","death")

fit.weib <- idm(formula02=Hist(time,event=death,entry=entry)~certif,
formula01=Hist(time=list(L,R),event=dementia)~certif,data=d
,irec=0,igraph=1,hazard="Splines",eps=c(1,1,1)) 
print(fit.weib) 
}
}

\keyword{methods}
