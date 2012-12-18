\name{illness}
\alias{illness}
\docType{data}
\title{
Data set for illness-death models
}
\description{
A simulated data frame for illness-death models composed of left-truncated, right-censored and interval-censored data.
}
\usage{data(illness)}
\format{
  A data frame with 3000 observations on the following 9 variables.
  \describe{
    \item{\code{indma}}{disease status, 0=non-diseased, 1=diseased}
    \item{\code{inddc}}{death status, 0=alive, 1=dead}
    \item{\code{t0}}{left-truncation time}
    \item{\code{t1}}{for diseased subjects: left endpoint of censoring interval; for non-diseased subjects: right censoring time for the 0 --> 1 transition}
    \item{\code{t2}}{for diseased subjects: right endpoint of censoring interval; for non-diseased subjects: right censoring time for the disease event}
    \item{\code{t3}}{for dead subjects: age at death; for alive subject: censoring time for the death event}
    \item{\code{var1}}{covariate 1}
    \item{\code{var2}}{covariate 2}
    \item{\code{var3}}{covariate 3}
  }
}

\examples{
data(illness)
}
\keyword{datasets}
