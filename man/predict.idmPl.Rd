\name{predict.idmPl}
\alias{predict.idmPl}

\title{
Predictions for an illness-death model using the penalized likelihood approach.
}
\description{
  Predict transition probabilities and cumulative probabilities from an object of class \code{idmPl}. Confidence intervals are calculated. 
}
\usage{
 \method{predict}{idmPl}(object, s, t, Z01, Z02, Z12, nsim=2000, CI=TRUE, ...)
}
\arguments{
  \item{object}{an \code{idmPl} class objects returned by a call to the \code{\link{idm}} function with \code{hazard}="Splines".}
  \item{s}{time at prediction.}
  \item{t}{time for prediction.}
  \item{Z01}{vector for the values of the covariates on the transition 0 --> 1 (in the same order as the covariates within the call. The default values are all 0.}  
  \item{Z02}{vector for the values of the covariates on the transition 0 --> 2 (in the same order as the covariates within the call. The default values are all 0.} 
  \item{Z12}{vector for the values of the covariates on the transition 1 --> 2 (in the same order as the covariates within the call. The default values are all 0.} 
  \item{nsim}{number of simulations for the confidence intervals calculations. The default is 2000.}
  \item{CI}{boolean: with (\code{TRUE}) or without (\code{WRONG}) confidence intervals for the predicted values. The default is \code{TRUE}.}
  \item{\dots}{other parameters.} 
}


\value{
     a list containing the following predictions with confidence bands:
   \item{p00}{the transition probability \eqn{p_{00}}.}
   \item{p01}{the transition probability \eqn{p_{01}}.}
   \item{p11}{the transition probability \eqn{p_{11}}.}
   \item{p12}{the transition probability \eqn{p_{12}}.}
   \item{p02_0}{the probability of direct transition from state 0 to state 2.}
   \item{p02_1}{the probability of transition from state 0 to state 2 via state 1.}
   \item{p02}{transition probability \eqn{p_{02}}. Note that \code{p02}=\code{p_02_0}+\code{p02_1}.}
   \item{F01}{the lifetime risk of disease. \code{F01}=\code{p01}+\code{p02_1}.}
   \item{F0.}{the probability of exit from state 0. \code{F0.}=\code{p02_0}+\code{p01}+\code{p02_1}.}
}

\seealso{
\code{\link{idm}}
}

\author{
J. Pierre
}

\examples{
\dontrun{
data(Dem15_h_3)
d <- Dem15_h_3
names(d) <- c("dementia","mort","entry","L","R","time","cep")

fit <- idm(formula02=Hist(time,event=mort,entry=entry)~cep,
formula01=Hist(time=list(L,R),event=dementia)~cep,data=d,hazard="Splines") 

pred <- predict(fit,s=70,t=80,Z01=c(1),Z02=c(1),Z12=c(1))
}
}

\keyword{methods}

