
\name{plot.survWeib}
\alias{plot.survWeib}


\title{Plot Method for survival using a parametrical approach.}

\description{
  Plot estimated survival hazard function, more generally of a class `survWeib' object. Confidence bands are allowed. 
}

\usage{

\method{plot}{survWeib}(x, ...)

}

\arguments{
  \item{x}{a fitted survival model using a parametrical weibull function (output from calling \code{surv}}
  \item{\dots}{ Other graphical parameters like those in \code{\link{plot.surv}}}

}

\value{
     Print a plot of hazard function and survival function of a class \code{survWeib} object
}

\seealso{\code{\link{print.survWeib}},\code{\link{summary.survWeib}},\code{\link{surv}},
}


\examples{


\dontrun{
/*** Weibull survival model ***/
data(test)
fit.su <- surv(Hist(time=list(l,r),id)~cov,data=test) 

# no confidence bands
plot(fit.su,conf.bands=FALSE)

# confidence bands
plot(fit.su)
 
} 

}

\keyword{methods}
