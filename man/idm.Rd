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
<<<<<<< HEAD
idm(formula01,
                formula02,
                formula12,
                data,
                maxiter=200,
                eps=c(5,5,3),
                n.knots=c(7,7,7),
                knots="equidistant",
                CV=FALSE,
                kappa=c(1000000,500000,20000),
                intensities="Weib",
		conf.int=TRUE,
                print.iter=FALSE,
                subset=NULL,
                na.action = na.fail)
}
\arguments{
  \item{formula01}{
A formula specifying a regression model for the \code{0 --> 1}
transition from the initial state to the transient state of the
illness-death model. 
The right hand side of the formula specifies the covariate terms, and the left hand side  must be an event
history object as returned by the function \code{Hist}. 
}
\item{formula02}{
A formula specifying a regression model for the \code{0 --> 2} transition
from the initial state to the absorbing state.
The left hand side must be equal to the left hand side of
  \code{formula01}. If missing it is set to \code{formula01}.
}
\item{formula12}{
  A formula specifying a regression model for the \code{1 --> 2} transition
  from the transient state to the absorbing state.
  operator is not required. If missing it is set to \code{formula01}.
}
  \item{data}{
A data frame in which to interpret the variables of \code{formula01}, \code{formula02} and \code{formula12}.
}
  \item{maxiter}{
Maximum number of iterations. The default is 200. 
}
  \item{eps}{
 A vector of 3 integers >0 used to define the power of three convergence criteria: 1. for the
  regression parameters, 2. for the likelihood, 3. for the second
  derivatives. The default is \code{c(5,5,3)} which is translated into 
  convergence if the respective values change less then
  \eqn{10^{-5}} (for regression parameters and likelihood) and
  \eqn{10^{-3}} for the second derivatives between two iterations.
}
\item{n.knots}{
  For \code{intensities="Splines"} only, a vector of length 3
  specifing the number of knots, one for each transition, for the M-splines estimate
  of the baseline intensities in the order
  \code{0 --> 1}, \code{0 --> 2}, \code{1 --> 2}.
  The default is c(7,7,7). 
}
\item{knots}{
List of length 3 containing the placements (timepoints) of the knots for
  the M-spline of the three transitions.  
}
  \item{CV}{
    Binary variable equals to 1 when search (by approximated cross validation) of the smoothing parameters kappa and 0 otherwise. Argument for the penalized likelihood approach. The default is 0.
}
  \item{kappa}{
a vector of length 3.  If CV=FALSE, smoothing parameters for the transition 0 --> 1, 0 --> 2 and 1 --> 2. If CV=TRUE, initial values of the smoothing parameters for the cross validation search. Argument for the penalized likelihood approach.  
}
  \item{intensities}{
type of estimation method: "Splines" for a penalized likelihood approach with approximation of the transition intensities by M-splines, "Weib" for a parametric approach with a Weibull distribution on the transition intensities. Default is "Weib".
}
  \item{conf.int}{
Boolean parameter. Equals to \code{TRUE} to calculate pointwise confidence intervals for the transition intensities curves, \code{FALSE} otherwise. Default is \code{TRUE}.
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
  \item{intensity01}{matched values of the intensities for transition 0 --> 1.} 
  \item{lowerIntensity01}{lower confidence intervals for the values of the intensities for transition 0 --> 1.}
  \item{upperIntensity01}{upper confidence intervals for the values of the intensities for transition 0 --> 1.}
  \item{intensity02}{matched values of the intensities for transition 0 --> 2.}
  \item{lowerIntensity02}{lower confidence intervals for the values of the intensities for transition 0 --> 2.} 
  \item{upperIntensity02}{upper confidence intervals for the values of the intensities for transition 0 --> 2.}
  \item{intensity12}{matched values of the intensities for transition 1 --> 2.}
  \item{lowerIntensity12}{lower confidence intervals for the values of the intensities for transition 1 --> 2.}
  \item{upperIntensity12}{upper confidence intervals for the values of the intensities for transition 1 --> 2.} 
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
  \item{CV}{a binary variable equals to 1 when search of the smoothing parameters \link{kappa} by approximated cross-validation, 1 otherwise. The default is 0.}    
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
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>

}

\seealso{
\code{\link{print.idmWeib}},
\code{\link{print.idmSplines}},
\code{\link{summary.idmWeib}},
\code{\link{summary.idmSplines}}
}

\examples{
\dontrun{

data(Paq1000)

# Illness-death model with certif on the 3 transitions
# Weibull parametrization and likelihood maximization
fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
formula01=Hist(time=list(l,r),event=dementia)~certif,data=Paq1000) 

fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
		formula01=Hist(time=list(l,r),event=dementia)~certif,
		data=Paq1000)

# Illness-death model with certif on transitions 01 and 02
# Splines parametrization and penalized likelihood maximization
fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=e)~certif,
		formula01=Hist(time=list(l,r),event=dementia)~certif,
                formula12=~1,
                intensities="Splines",
		data=Paq1000)

## to print
fit.weib

## to summary
summary(fit.splines)
}
}

\keyword{ilness-death}

