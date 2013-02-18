\name{Paq1000}
\alias{Paq1000}
\docType{data}
\title{
Paquid data set
}
\description{
Paquid data set composed of 1000 subjects selected randomly from the Paquid data set of 3675 subjects.
}
\usage{data(Paq1000)}
\format{
  A data frame with 1000 rows and the following 8 columns.
  \describe{
    \item{\code{dementia}}{dementia status, 0=non-demented, 1=demented}
    \item{\code{death}}{death status, 0=alive, 1=dead}
    \item{\code{t0}}{age at entry in the study}
    \item{\code{t1}}{for demented subjects: age at the visit before the diagnostic visit; for non-demented subjects: age at the last visit (censoring age)}
    \item{\code{t2}}{for demented subjects: age at the diagnostic visit; for non-demented subjects: age at the last visit (censoring age)}
    \item{\code{t3}}{for dead subjects: age at death; for alive subject: age at the latest news}
    \item{\code{certif}}{primary school certificate}
    \item{\code{gender}}{gender, 0=female, 1=male}
  }
}

\examples{
data(Paq1000)
head(1000)
}
\keyword{datasets}
