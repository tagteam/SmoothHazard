<<<<<<< HEAD
\name{plot.idm}
\alias{plot.idm}
\title{Plot method for an illness-death model}
\usage{
\method{plot}{idm}(x, conf.int = FALSE, citype = "shadow", add = FALSE,
  axes = TRUE, col, lwd, lty, xlim, ylim, xlab, ylab, legend = TRUE,
  transition = c("01", "02", "12"), ...)
}
\arguments{
  \item{x}{a \code{idmWeib} class object (output from
  calling \code{idm} with the (default) option
  \code{intensities}="Weib".}

  \item{conf.int}{If TRUE show confidence limits}

  \item{citype}{Type of confidence limits, can be "shadow"
  or "bars"}

  \item{add}{If TRUE add to existing plot}

  \item{axes}{If TRUE axes are drawn}

  \item{col}{Color of the lines}

  \item{lwd}{Width of the lines}

  \item{lty}{Type of the lines}

  \item{xlim}{Limits for x-axis}

  \item{ylim}{Limits for y-axis}

  \item{xlab}{Label for x-axis}

  \item{ylab}{Label for y-axis}

  \item{legend}{If TRUE a legend is drawn, which can be
  further controlled via \code{\link{SmartControl}}.}

  \item{transition}{Choose one of the transition
  intensities: \code{c("01","02","12")}.}

  \item{...}{Passed to \code{\link{SmartControl}}}
}
\value{
Print a plot of the baseline transition intensities of an
illness-death model estimated using a Weibull approach.
}
\description{
Plot estimated baseline transition intensities from an
object of class \code{idm} optionally with confidence
limits.
}
\examples{
library(lava)
m <- idmModel(scale.lifetime=1/10,scale.illtime=1/8)
distribution(m,"X") <- binomial.lvm()
regression(m,latent.lifetime~X) <- 0.7
set.seed(30)
d <- sim(m,100)
fit.weib <- idm(formula02=Hist(observed.lifetime,event=seen.exit)~1,
formula01=Hist(time=list(L,R),event=seen.ill)~1,data=d,conf.int=FALSE)
plot(fit.weib)

\dontrun{
## FIXME: the limits for the 01 transition are a bit wide!?
## with bootstrap confidence limits
fit.weib <- idm(formula02=Hist(observed.lifetime,event=seen.exit)~1,
formula01=Hist(time=list(L,R),event=seen.ill)~1,data=d,conf.int=TRUE)
plot(fit.weib)
}
}
\seealso{
\code{\link{print.idm}},\code{\link{summary.idm}},\code{\link{idm}},

\code{\link{idm}}
}
\keyword{methods}

=======
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

\method{plot}{idm}(x, conf.int=FALSE, ylim,
 pos.legend = "topleft", main, transition, ...)

}
\arguments{
  \item{x}{a \code{idmWeib} class object or a \code{idmPl} class object (output from calling \code{\link{idm}}
    function).}
  \item{conf.int}{logical value. Determines whether pointwise confidence
    intervals will be plotted. The default is FALSE.}
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

# no pointwise confidence intervals
plot(fit.weib)

# pointwise confidence intervals for transition 0 --> 1
plot(fit.weib,conf.int=TRUE,transition="01")
}
}

\keyword{methods}

>>>>>>> saves before pulling Celia's changes
