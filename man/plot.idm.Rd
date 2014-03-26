\name{plot.idm}
\alias{plot.idm}

\title{
Plot method for an illness-death model.
}
\description{
  Plot estimated baseline transition intensities from an object of class
  \code{idmWeib} or \code{idmPl}. Pointwise confidence intervals are available. 
}
\usage{
\method{plot}{idm}(x,
		    conf.int=FALSE,
		    citype="shadow",
		    add=FALSE,
		    axes=TRUE,
		    col,
		    lwd,
		    lty,
		    xlim,
		    ylim,
		    xlab,
		    ylab,
		    legend=TRUE,
		    transition=c("01","02","12"),
		    ...)
}
\arguments{
  \item{x}{a \code{idmWeib} class object or a \code{idmSplines} class object (output from calling \code{\link{idm}}
    function).}
  \item{conf.int}{logical value. Determines whether pointwise confidence
    intervals will be plotted. The default is FALSE.}
  \item{citype}{Type of pointwise confidence intervals, either "shadow"
    or "lines".}
  \item{add}{If TRUE no new plot is generated.}
  \item{axes}{If FALSE no axes are drawn.}
  \item{col}{Colour of the lines}
  \item{lwd}{Thickness of the lines}
  \item{lty}{lty of the lines}
  \item{xlim}{Limits of the x-axis}
  \item{ylim}{vector containing the y limites (min, max) of the
    plot. The default is the min and the max of all the values of the
    three transition intensities}  
  \item{xlab}{Label for the x-axis}
  \item{ylab}{Label for the y-axis}
  \item{legend}{Whether to show the legend.}
  \item{ylim}{vector containing the y limites (min, max) of the plot. The default is the min and the max of all the values of the three transition intensities}  
 \item{pos.legend}{The location of the legend can be specified by setting this argument to a single keyword from the  list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and
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
  Thomas Alexander. Gerds
}

\examples{
\dontrun{
data(Paq1000)
fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
		formula01=Hist(time=list(l,r),event=dementia)~certif,
		data=Paq1000) 

# no pointwise confidence intervals
plot(fit.weib)

# pointwise confidence intervals for transition 0 --> 1
plot(fit.weib,conf.int=TRUE,transition="01")
}
}

\keyword{methods}

