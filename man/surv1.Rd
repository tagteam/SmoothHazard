\name{surv1}
\alias{surv1}
\docType{data}
\title{
Data set for survival models: left-truncated, right-censored and interval-censored data
}
\description{
 A simulated data frame for survival models composed of left-truncated, right-censored and interval-censored data.
}
\usage{data(surv1)}
\format{
  A data frame with 3000 observations on the following 7 variables.
  \describe{
    \item{\code{ind}}{disease status}
    \item{\code{t0}}{left-truncation time}
    \item{\code{t1}}{for diseased subjects: left endpoint of censoring interval; for non-diseased subjects: right censoring time}
    \item{\code{t2}}{for diseased subjects: right endpoint of censoring interval; for non-diseased subjects: right censoring time for the disease event}
    \item{\code{x1}}{covariate 1}
    \item{\code{x2}}{covariate 2}
    \item{\code{x3}}{covariate 3}
  }
}
\examples{
data(surv1)
}
\keyword{datasets}
