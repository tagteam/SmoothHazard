
\name{plot.survPl}
\alias{plot.survPl}


\title{Plot method for a survival model using a penalized likelihood approach.}

\description{
  Plot estimated baseline survival function from an object of class \code{survPl}. Confidence bands are allowed. 
}

\usage{

\method{plot}{survPl}(x, ...)

}

\arguments{
  \item{x}{a fitted survival model using a penalized likelihood approach (output from calling \code{surv} with the option \code{hazard}="Spline".}
  \item{\dots}{other graphical parameters like those in \code{\link{plot.surv}}}
}

\value{
     Print a plot of hazard function and survival function of a class \code{survPl} object
}

\seealso{\code{\link{print.survPl}},\code{\link{summary.survPl}},\code{\link{surv}},
}


\examples{


\dontrun{
/*** a penalized survival model ***/
data(test)
fit.su <- surv(Hist(time=list(l,r),id)~cov,data=test,hazard="Splines") 

# no confidence bands
plot(fit.su,conf.bands=FALSE)

# confidence bands
plot(fit.su)
 
} 

}

\keyword{methods}
