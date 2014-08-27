\name{plot.idmSplines}
\alias{plot.idmSplines}

\title{Plot method for an illness-death model using a penalized likelihood approach.}

\description{
  Plot estimated baseline transition intensities from a
  \code{idmSplines} class object. Pointwise confidence intervals are available. 
}

\usage{

\method{plot}{idmSplines}(x, ...)

}

\arguments{
  \item{x}{a \code{idmSplines} class object (output from calling \code{idm} with the option \code{intensities}="Splines".}
  \item{\dots}{other graphical parameters like those in \code{\link{plot.idm}}}

}

\value{
     Print a plot of the baseline transition intensities of an illness-death model estimated using a penalized log-likelihood approach.
}

\seealso{\code{\link{print.idmSplines}},\code{\link{summary.idmSplines}},\code{\link{idm}},
}


\examples{


\dontrun{
# illness deaths model  ***/
data(Paq1000)

fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=e)~certif,
		formula01=Hist(time=list(l,r),event=dementia)~certif,
                formula12=~1,
                intensities="Splines",
		data=Paq1000)
plot(fit.splines)
 
} 

}

\keyword{methods}
