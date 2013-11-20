
\name{plot.shrSplines}
\alias{plot.shrSplines}


\title{Plot method for a survival model using a penalized likelihood approach.}

\description{
  Plot estimated baseline survival function from an object of class
  \code{shrSplines}. Pointwise confidence limits are available. 
}

\usage{

\method{plot}{shrSplines}(x, ...)

}

\arguments{
  \item{x}{a fitted survival model using a penalized likelihood approach (output from calling \code{shr} with the option \code{hazard}="Splines".}
  \item{\dots}{other graphical parameters like those in \code{\link{plot.shr}}}
}

\value{
     Print a plot of hazard function and survival function of a class \code{shrSplines} object
}

\seealso{\code{\link{print.shrSplines}},\code{\link{summary.shrSplines}},\code{\link{shr}},
}


\examples{


\dontrun{
/*** a penalized survival model ***/
data(test)
fit.su <- shr(Hist(time=list(l,r),id)~cov,data=test,hazard="Splines") 

# no pointwise confidence limits
plot(fit.su,conf.int=FALSE)

# pointwise confidence limits
plot(fit.su)
 
} 

}

\keyword{methods}
