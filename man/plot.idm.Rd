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
 pos.legend = "topleft", main, transition, ...)

}
\arguments{
  \item{x}{a \code{idmWeib} class object or a \code{idmPl} class object (output from calling \code{\link{idm}}
    function).}
  \item{conf.int}{logical value. Determines whether confidence bands will be plotted. The default is FALSE.}
  \item{ylim}{vector containing the y limites (min, max) of the plot. The default is the min and the max of all the values of the three transition intensities}  
 \item{pos.legend}{The location of the legend can be specified by setting this argument to a single keyword from the 
 list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and
     "center". The default is "topleft"}
  \item{main}{title of plot}
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
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>

}

\examples{
\dontrun{
# Weibull illness deaths model 
data(Dem15_h_3)
d <- Dem15_h_3
names(d) <- c("dementia","mort","entry","L","R","time","cep")

fit.weib <- idm(formula02=Hist(time,event=mort,entry=entry)~cep,
formula01=Hist(time=list(L,R),event=dementia)~cep,data=d) 

# no confidence bands
plot(fit.weib)

# confidence bands for transition 0 --> 1
plot(fit.weib,conf.int=TRUE,transition="01")
}
}

\keyword{methods}

