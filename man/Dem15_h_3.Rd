\name{Dem15_h_3}
\alias{Dem15_h_3}
\docType{data}
\title{
Paquid data set
}
\description{
Paquid hommes suivi 15 ans.
}
\usage{data(Dem15_h_3)}
\format{
  A data frame with 1540 rows and the following 7 columns.
  \describe{
    \item{\code{dementia}}{dementia status, 0=non-demented, 1=demented}
    \item{\code{mort}}{death status, 0=alive, 1=dead}
    \item{\code{t0}}{age at entry in the study}
    \item{\code{t1}}{for demented subjects: age at the visit before the diagnostic visit; for non-demented subjects: age at the last visit (censoring age)}
    \item{\code{t2}}{for demented subjects: age at the diagnostic visit; for non-demented subjects: age at the last visit (censoring age)}
    \item{\code{t3}}{for dead subjects: age at death; for alive subject: age at the latest news}
    \item{\code{cep}}{primary school certificate}
  }
}

\examples{
data(Dem15_h_3)
head(Dem15_h_3)
}
\keyword{datasets}
