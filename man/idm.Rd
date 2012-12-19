\name{illness-death model}
\alias{idm}

\title{
 Fit an illness-death model
}
\description{
Fit an illness-death model using either a semi-parametric approach (penalized likelihood with an approximation of the transition intensity functions by linear combination of M-splines) or a parametric approach (specifying Weibull distributions on the transition intensities).
Left-truncated, right-censored, and interval-censored data are allowed.
State 0 corresponds to the initial state, state 1 to the transient one, state 2 to the absorbant one. The allowed transitions are: 0 --> 1, 0 --> 2 and 1 --> 2. 
}
\usage{
idm(formula01, formula02, formula12, data, maxiter=200, eps=c(5,5,3), 
n.knots=c(7,7,7), irec=0,kappa0 = c(1000000,500000,20000), igraph=1, 
hazard="Weib", print.iter=FALSE, subset=NULL, na.action=na.omit)
}

\arguments{
  \item{formula01}{
a formula object for transition '0 --> 1' (from the initial state to the transient state), with the response on the left of a \eqn{\texttildelow}
operator, and the terms on the right. The response must be a
survival object or Hist object as returned by the 'Surv' or 'Hist' function.
}
  \item{formula02}{
a formula object for transition '0 --> 2' (from the initial state to the absorbant state), with the response on the left of a \eqn{\texttildelow}
operator, and the terms on the right. The response must be a
survival object or Hist object as returned by the 'Surv' or 'Hist' function.
}
  \item{formula12}{
a formula object for transition '1 --> 2' (from the transient state to the absorbant state). The response on the left of the \eqn{\texttildelow}
operator is not required. If missing formula12, variables of transition '1--->2' are the same as the ones of transition '0--->1'.
}
  \item{data}{
a data frame in which to interpret the variables named in \code{formula01}, \code{formula02} and \code{formula12}.
}
  \item{maxiter}{
maximum number of iterations. The default is 200. 
}
  \item{eps}{
a vector of length 3 for the convergence criteria (criterion for parameters, criterion for likelihood, criterion for second derivatives). The default is 'c(5,5,3)' and corresponds to criteria equals to \eqn{10^{-5}}, \eqn{10^{-5}} and \eqn{10^{-3}}.
}
  \item{n.knots}{
a vector of length 3 for the number of knots for the splines to use to approximate the intensities for transitions 0 --> 1, 0 --> 2 and 1 --> 2 respectively. Argument for the penalized likelihood approach.
  The default is c(7,7,7). 
}
  \item{irec}{
binary variable equals to 1 when search (by approximated cross validation) of the smoothing parameters kappa and 0 otherwise. Argument for the penalized likelihood approach. The default is 0.
}
  \item{kappa0}{
a vector of length 3.  If irec=0, smoothing parameters for the transition 0 --> 1, 0 --> 2 and 1 --> 2. If irec=1, initial values of the smoothing parameters for the cross validation search. Argument for the penalized likelihood approach.  
}
  \item{igraph}{
a binary variable equals to 1 when you want to output graphs of the transition intensities, 0 otherwise. The default is 1.
}
  \item{hazard}{
type of estimation method: "Splines" for a penalized likelihood approach with approximation of the transition intensities by M-splines, "Weib" for a parametric approach with a Weibull distribution on the transition intensities. Default is "Weib".
}
  \item{print.iter}{
boolean parameter. Equals to \code{TRUE} to print the likelihood during the iteration process, \code{FALSE} otherwise. Default is \code{FALSE}. This option is not running on Windows.
}
  \item{subset}{
expression indicating the subset of the rows of data to be used in the fit. All observaation are included by default.
}
  \item{na.action}{
how NAs are  treated. The default is first, any na.action attribute of data, second a na.action setting of options, and third 'na.fail' if that is unset. The 'factory-fresh' default is na.omit. Another possible value is NULL.
}
}
\details{
  The estimated parameters are obtained using the robust Marquardt algorithm (Marquardt, 1963) which is a combination
  between a Newton-Raphson algorithm and a steepest descent algorithm.
}
\value{

  \item{call}{the call that produced the result.}
  \item{coef}{regression parameters.}
  \item{loglik}{vector containing the log-likelihood without and with covariate.}
  \item{cv}{vector containing the convergence criteria.}
  \item{niter}{number of iterations.} 
  \item{converged}{integer equal to 1 when the model converged, 2, 3 or 4 otherwise.}
  \item{modelPar}{Weibull parameters.}
  \item{N}{number of subjects.}
  \item{events1}{number of events 0 --> 1.}
  \item{events2}{number of events 0 --> 2 or 0 --> 1 --> 2.}
  \item{NC}{vector containing the number of covariates on transitions 0 --> 1, 0 --> 2, 1 --> 2.} 
  \item{responseTrans}{model response for the 0 --> 1 transition. \code{Hist} or \code{Surv} object.}
  \item{responseAbs}{model response for the 0 --> 2 transition. \code{Hist} or \code{Surv} object.}
  \item{time}{times for which transition intensities have been evaluated for plotting. Vector in the Weibull approach. Matrix in the penalized likelihhod approach for which the colums corresponds to the transitions 0 --> 1, 1 --> 2, 0 --> 2.}
  \item{hazard01}{matched values of the intensities for transition 0 --> 1.} 
  \item{lowerHazard01}{lower confidence intervals for the values of the intensities for transition 0 --> 1.}
  \item{upperHazard01}{upper confidence intervals for the values of the intensities for transition 0 --> 1.}
  \item{hazard02}{matched values of the intensities for transition 0 --> 2.}
  \item{lowerHazard02}{lower confidence intervals for the values of the intensities for transition 0 --> 2.} 
  \item{upperHazard02}{upper confidence intervals for the values of the intensities for transition 0 --> 2.}
  \item{hazard12}{matched values of the intensities for transition 1 --> 2.}
  \item{lowerHazard12}{lower confidence intervals for the values of the intensities for transition 1 --> 2.}
  \item{upperHazard12}{upper confidence intervals for the values of the intensities for transition 1 --> 2.} 
  \item{RR}{vector of relative risks.}
  \item{V}{variance-covariance matrix.}
  \item{se}{standart errors of the regression parameters.}
  \item{Xnames01}{names of covariates on 0 --> 1.}
  \item{Xnames02}{names of covariates on 0 --> 2.}
  \item{Xnames12}{names of covariates on 1 --> 2.}
  \item{knots01}{knots to approximate by M-splines the intensity of the 0 --> 1 transition.} 
  \item{knots02}{knots to approximate by M-splines the intensity of the 0 --> 2 transition.}
  \item{knots12}{knots to approximate by M-splines the intensity of the 1 --> 2 transition.}
  \item{nknots01}{number of knots on transition 0 --> 1.} 
  \item{nknots02}{number of knots on transition 0 --> 2.}  
  \item{nknots12}{number of knots on transition 1 --> 2.}    
  \item{theta01}{square root of splines coefficients for transition 0 --> 1.}
  \item{theta02}{square root of splines coefficients for transition 0 --> 2.} 
  \item{theta12}{square root of splines coefficients for transition 1 --> 2.}   
  \item{irec}{a binary variable equals to 1 when search of the smoothing parameters \link{kappa} by approximated cross-validation, 1 otherwise. The default is 0.}    
  \item{igraph}{a binary variable equals 1 when you need to output transition intensity to plot, 0 otherwise. The default is 1.} 
  \item{kappa}{vector containing the smoothing parameters for transition 0 --> 1, 0 --> 2, 1 --> 2 used to estimate the model by the penalized likelihood approach.} 
  \item{CVcrit}{cross validation criteria.}    
  \item{DoF}{degrees of freedom of the model.} 
  \item{na.action}{observations deleted if missing values.} 
}

\references{
D. Marquardt (1963). An algorithm for least-squares estimation of nonlinear parameters. 
\emph{SIAM Journal of Applied Mathematics}, 431-441. 
}

\author{
J. Pierre
}

\seealso{
\code{\link{print.idmWeib}},
\code{\link{print.idmPl}},
\code{\link{summary.idmWeib}},
\code{\link{summary.idmPl}}
}

\examples{
# Weibull illness deaths model 
data(Dem15_h_3)
d <- Dem15_h_3
names(d) <- c("dementia","mort","entry","L","R","time","cep")

fit.weib <- idm(formula02=Hist(time,event=mort,entry=entry)~cep,
		formula01=Hist(time=list(L,R),event=dementia)~cep,
		data=d)


fit.weib <- idm(formula02=Hist(time,event=mort,entry=entry)~cep,
formula01=Hist(time=list(L,R),event=dementia)~cep,data=d) 

## to print
fit.weib

## to summary
summary(fit.weib)

}

\keyword{ilness-death}

