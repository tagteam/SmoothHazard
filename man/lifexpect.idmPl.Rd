\name{lifexpect.idmPl}
\alias{lifexpect.idmPl}

\title{
Predictions of life expectancies from an illness-death model using a penalized likelihood approach.
}

\description{
  Predict life expectancies from an object of class \code{idmPl}. Life expectancies are calculated at time \code{s} for a subject who has the covariates values Z01, Z02, Z12. Confidence intervals are calculated. 
}

\usage{lifexpect.idmPl(object, s, Z01, Z02, Z12, nsim=1000, CI=TRUE, 
...)
}

\arguments{
  \item{object}{an \code{idmPl} class objects returned by a call to the \code{\link{idm}} function with \code{hazard}="Splines".}
  \item{s}{time at prediction.}
  \item{Z01}{vector for the values of the covariates on the transition 0 --> 1 (in the same order as the covariates within the call. The default values are all 0.}  
  \item{Z02}{vector for the values of the covariates on the transition 0 --> 2 (in the same order as the covariates within the call. The default values are all 0.}
  \item{Z12}{vector for the values of the covariates on the transition 1 --> 2 (in the same order as the covariates within the call. The default values are all 0.}    
  \item{nsim}{number of simulations for the confidence intervals calculations. The default is 1000.}
  \item{CI}{boolean: with (\code{TRUE}) or without (\code{WRONG}) confidence intervals for the life expectancies. The default is \code{TRUE}.}
  \item{\dots}{others arguments.}
}

\value{
   a list containing:
   \item{HLE}{healthy life expectancy and confidence interval.}
   \item{LE0}{life expectancy of a non-diseased subject and confidence interval.} 
   \item{LE1}{life expectancy of a diseased subject and confidence interval.}
}

\seealso{
\code{\link{idm}}
}

\author{
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>

}

\examples{
\dontrun{
data(Paq1000)

fit <- idm(formula02=Hist(time,event=death,entry=entry)~certif,
formula01=Hist(time=list(L,R),event=dementia)~certif,data=d,hazard="Splines") 

LE <- lifexpect.idmPl(fit,s=95,Z01=c(1),Z02=c(1),Z12=c(1))
}
}

\keyword{methods}

