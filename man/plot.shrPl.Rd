
\name{plot.shrPl}
\alias{plot.shrPl}


\title{Plot method for a survival model using a penalized likelihood approach.}

\description{
  Plot estimated baseline survival function from an object of class \code{shrPl}. Confidence bands are allowed. 
}

\usage{

\method{plot}{shrPl}(x, ...)

}

\arguments{
  \item{x}{a fitted survival model using a penalized likelihood approach (output from calling \code{shr} with the option \code{hazard}="Spline".}
  \item{\dots}{other graphical parameters like those in \code{\link{plot.shr}}}
}

\value{
     Print a plot of hazard function and survival function of a class \code{shrPl} object
}

\seealso{\code{\link{print.shrPl}},\code{\link{summary.shrPl}},\code{\link{shr}},
}


\examples{


\dontrun{
/*** a penalized survival model ***/
data(test)
fit.su <- shr(Hist(time=list(l,r),id)~cov,data=test,hazard="Splines") 

# no confidence bands
plot(fit.su,conf.bands=FALSE)

# confidence bands
plot(fit.su)
 
} 

}

\keyword{methods}
