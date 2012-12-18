\name{test}
\alias{test}
\docType{data}
\title{
Data set for survival models: right-censored and interval-censored data.
}
\description{
A simulated data frame for survival models composed of right-censored and interval-censored data.
}
\usage{data(test)}
\format{
  A data frame with 936 observations on the following 4 variables.
  \describe{
    \item{\code{l}}{for diseased subjects: left endpoint of censoring interval; for non-diseased subjects: right censoring time}
    \item{\code{r}}{for diseased subjects: right endpoint of censoring interval; for non-diseased subjects: right censoring time for the disease event}
    \item{\code{id}}{disease status}
    \item{\code{cov}}{covariate}
  }
}

\examples{
data(test)
head(test)
}
\keyword{datasets}
