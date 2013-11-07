
\name{plot.idmWeib}
\alias{plot.idmWeib}


\title{Plot method for an illness-death model using a Weibull approach.}

\description{
  Plot estimated baseline transition intensities from an object of class \code{idmWeib}. Confidence bands are allowed. 
}

\usage{

\method{plot}{idmWeib}(x, ...)

}

\arguments{
  \item{x}{a \code{idmWeib} class object (output from calling \code{idm} with the (default) option \code{intensities}="Weib".}
  \item{\dots}{other graphical parameters like those in \code{\link{plot.idm}}.}

}

\value{
     Print a plot of the baseline transition intensities of an illness-death model estimated using a Weibull approach.
}

\seealso{\code{\link{print.idmWeib}},\code{\link{summary.idmWeib}},\code{\link{idm}},
}


\examples{


\dontrun{
data(Paq1000)

fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~1,
formula01=Hist(time=list(l,r),event=dementia)~1,data=Paq1000)

plot(fit.weib)
 
} 

}

\keyword{methods}
