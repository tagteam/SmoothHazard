
\name{plot.shrWeib}
\alias{plot.shrWeib}


\title{Plot Method for survival using a parametrical approach.}

\description{
  Plot estimated survival hazard function, more generally of a class
  `shrWeib' object. Pointwise confidence limits are available.
}

\usage{

\method{plot}{shrWeib}(x, ...)

}

\arguments{
  \item{x}{a fitted survival model using a parametrical weibull function (output from calling \code{shr}}
  \item{\dots}{ Other graphical parameters like those in \code{\link{plot.shr}}}

}

\value{
     Print a plot of hazard function and survival function of a class \code{shrWeib} object
}

\seealso{\code{\link{print.shrWeib}},\code{\link{summary.shrWeib}},\code{\link{shr}},
}


\examples{


\dontrun{
/*** Weibull survival model ***/
data(testdata)
fit.su <- shr(Hist(time=list(l,r),id)~cov,data=testdata) 

# no pointwise confidence limits
plot(fit.su,conf.int=FALSE)

# pointwise confidence limits
plot(fit.su)
 
} 

}

\keyword{methods}
