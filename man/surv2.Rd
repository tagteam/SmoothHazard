\name{surv2}
\alias{surv2}
\docType{data}
\title{
Data set for survival models: left-truncated and right-censored data.
}
\description{
A simulated data frame for survival models composed of left-truncated and right-censored data.
}
\usage{data(surv2)}
\format{
  A data frame with 3000 observations on the following 6 variables.
  \describe{
    \item{\code{ind}}{disease status}
    \item{\code{t0}}{left-truncation time}
    \item{\code{t}}{for diseased subjects: exact event time; for non-diseased subjects: right censoring time}
    \item{\code{x1}}{covariate 1}
    \item{\code{x2}}{covariate 2}
    \item{\code{x3}}{covariate 3}
  }
}
\examples{
data(surv2)
}
\keyword{datasets}
