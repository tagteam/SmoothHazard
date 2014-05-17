\name{lifexpect}
\alias{lifexpect}
\title{Predictions of life expectancies from an illness-death.}
\usage{
lifexpect(x, s, Z01, Z02, Z12, nsim = 1000, CI = TRUE, ...)
}
\arguments{
  \item{x}{an object as returned by a call to the
  \code{\link{idm}} function.}

  \item{s}{time at prediction.}

  \item{Z01}{vector for the values of the covariates on the
  transition 0 --> 1 (in the same order as the covariates
  within the call. The default values are all 0.}

  \item{Z02}{vector for the values of the covariates on the
  transition 0 --> 2 (in the same order as the covariates
  within the call. The default values are all 0.}

  \item{Z12}{vector for the values of the covariates on the
  transition 1 --> 2 (in the same order as the covariates
  within the call. The default values are all 0.}

  \item{nsim}{number of simulations for the confidence
  intervals calculations. The default is 1000.}

  \item{CI}{boolean: with (\code{TRUE}) or without
  (\code{WRONG}) confidence intervals for the life
  expectancies. The default is \code{TRUE}.}

  \item{\dots}{others arguments.}
}
\value{
a list containing: \item{life.in.0.expectancy}{life
expectancy in state 0 and confidence interval.}
\item{life.expectancy.nondis}{life expectancy of a
non-diseased subject and confidence interval.}
\item{life.expectancy.dis}{life expectancy of a diseased
subject and confidence interval.}
}
\description{
Predict life expectancies from an object of class
\code{idm}. Life expectancies are calculated at time
\code{s} for a subject who has the covariates values Z01,
Z02, Z12. Confidence intervals are calculated.
}
\examples{
\dontrun{
data(Paq1000)

fit <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
formula01=Hist(time=list(l,r),event=dementia)~certif,method="Splines",data=Paq1000)

pred <- lifexpect(fit,s=70,t=80,Z01=c(1),Z02=c(1),Z12=c(1))

}
}
\author{
C. Touraine
}
\seealso{
\code{\link{idm}}
}
\keyword{methods}

