\name{plot.idmPl}
\alias{plot.idmPl}

\title{Plot method for an illness-death model using a penalized likelihood approach.}

\description{
  Plot estimated baseline transition intensities from a \code{idmPl} class object. Confidence bands are allowed. 
}

\usage{

\method{plot}{idmPl}(x, ...)

}

\arguments{
  \item{x}{a \code{idmPl} class object (output from calling \code{idm} with the option \code{hazard}="Splines".}
  \item{\dots}{other graphical parameters like those in \code{\link{plot.idm}}}

}

\value{
     Print a plot of the baseline transition intensities of an illness-death model estimated using a penalized log-likelihood approach.
}

\seealso{\code{\link{print.idmPl}},\code{\link{summary.idmPl}},\code{\link{idm}},
}


\examples{


\dontrun{
# illness deaths model  ***/
data(Paq1000)
d <- Paq1000
names(d) <- c("dementia","mort","entry","L","R","time","certif","death")

fit.weib <- idm(formula02=Hist(time,event=death,entry=entry)~certif,
formula01=Hist(time=list(L,R),event=dementia)~certif,data=d
,irec=0,igraph=1,hazard="Splines",eps=c(1,1,1)) 

# no confidence bands
plot(fit.weib)
 
} 

}

\keyword{methods}
