\name{plot.surv}
\alias{plot.surv}

\title{
Plot method for a survival model.
}
\description{
  Plot estimated baseline survival function from an object of class \code{surv}. Confidence bands are allowed. 
}
\usage{

\method{plot}{surv}(x, type="surv", add = FALSE, newdata=NULL, cause=NULL, col, lty, lwd,
	ylim, xlim, xlab="Time" ,ylab, legend=TRUE, confint=TRUE, timeOrigin=0,
	axes=TRUE, percent=TRUE,...)
}
\arguments{
  \item{x}{a \code{survWeib} or a \code{survPl} class object (output from calling \code{\link{surv}}
    function).
}
  \item{type}{type of function to plot. The default is "surv".
} 
  \item{add}{boolean.
}
  \item{newdata}{newdata.
}
  \item{cause}{cause.
}
  \item{col}{col.
}
  \item{lty}{lty.
}
  \item{lwd}{lwd.
}
  \item{ylim}{ylim.
}
  \item{xlim}{xlim.
}
  \item{xlab}{xlab.
}
  \item{ylab}{ylab.
}
  \item{legend}{legend.
}
  \item{confint}{confint.
}
  \item{timeOrigin}{timeOrigin.
}
  \item{axes}{axes.
}
  \item{percent}{percent.
}
  \item{\dots}{other graphical parameters.}
}


\value{
     Print a plot of a suvival model.
}

\seealso{
\code{\link{plot.surv}}
}

\author{
J. Pierre
}

\examples{
\dontrun{
# Weibull survival model 
data(test)
fit.su <- surv(Hist(time=list(l,r),id)~cov,data=test) 

# confidence bands
plot(fit.su)

# no confidence bands
plot(fit.su,confint=FALSE)

}
}

\keyword{methods}

