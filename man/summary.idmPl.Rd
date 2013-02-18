\name{summary.idmPl}
\alias{summary.idmPl}

\title{
Summary of a fitted illness-death model using a penalized likelihood approach
}

\usage{
\method{summary}{idmPl}(object, conf.int = 0.95, digits = 4,
 pvalDigits = 4, eps = 1e-04, ...)
}
\description{Print a short summary of a fitted illness-death model using the penalized likelihood approach.}
\arguments{
  \item{object}{
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
 \code{\link{idm}}, 
 \code{\link{print.idmPl}}, 
 \code{\link{plot.idm}} 
 \code{\link{plot.idmPl}}
}
  
\author{
J. Pierre
}

\examples{
\dontrun{
data(Paq1000)
d <- Paq1000
names(d) <- c("dementia","mort","entry","L","R","time","certif","death")

fit.weib <- idm(formula02=Hist(time,event=death,entry=entry)~certif,
formula01=Hist(time=list(L,R),event=dementia)~certif,data=d
,irec=0,igraph=1,hazard="Splines",eps=c(1,1,1)) 
summary(fit.weib) 
}
}

\keyword{methods}
