\name{plot.idm}
\alias{plot.idm}

\title{
Plot method for an illness-death model.
}
\description{
  Plot estimated baseline transition intensities from an object of class \code{idmWeib} or \code{idmPl}. Confidence bands are allowed. 
}
\usage{

\method{plot}{idm}(x, conf.int=FALSE, ylim,
transition, ...)

}
\arguments{
  \item{x}{a \code{idmWeib} class object or a \code{idmSplines} class object (output from calling \code{\link{idm}}
    function).}
  \item{conf.int}{logical value. Determines whether confidence bands will be plotted. The default is FALSE.}
  \item{ylim}{vector containing the y limites (min, max) of the plot. The default is the min and the max of all the values of the three transition intensities}  
  \item{transition}{a vector of the transitions to plotted on separate graphs. Possible value are "01","02" and "12". The default is a plot of the 3 intensities on the same graph.}  
  \item{\dots}{other graphical parameters.}
}


\value{
     Print a plot of the baseline transition intensities of an illness-death model.
}

\seealso{
\code{\link{idm}}
}

\author{
T. A. Gerds
}

\examples{
\dontrun{
data(Paq1000)
fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
formula01=Hist(time=list(l,r),event=dementia)~certif,data=Paq1000) 

# no confidence bands
plot(fit.weib)

# confidence bands for transition 0 --> 1
plot(fit.weib,conf.int=TRUE,transition="01")
}
}

\keyword{methods}

