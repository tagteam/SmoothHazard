
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
  \item{x}{a \code{idmWeib} class object (output from calling \code{idm} with the (default) option \code{hazard}="Weib".}
  \item{\dots}{other graphical parameters like those in \code{\link{plot.idm}}.}

}

\value{
     Print a plot of the baseline transition intensities of an illness-death model estimated using a Weibull approach.
}

\seealso{\code{\link{print.idmWeib}},\code{\link{summary.idmWeib}},\code{\link{idm}},
}


\examples{


\dontrun{
/*** Weibull illness deaths model ***/
data(Paq1000)
d <- Paq1000
names(d) <- c("dementia","mort","entry","L","R","time","certif","death")

fit.weib <- idm(formula02=Hist(time,event=death,entry=entry)~certif,
formula01=Hist(time=list(L,R),event=dementia)~certif,data=d) 

# no confidence bands
plot(fit.weib)
 
} 

}

\keyword{methods}
