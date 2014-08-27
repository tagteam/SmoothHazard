<<<<<<< HEAD
\name{shr}
\alias{shr}
\title{Fit a survival model}
\usage{
shr(formula, data, eps = c(5, 5, 3), nknots = 7, CV = FALSE,
  kappa = 10000, conf.int = TRUE, maxiter = 200, method = "Weib",
  print.iter = FALSE, na.action = na.omit)
}
\arguments{
  \item{formula}{a formula object with the response on the
  left of a \eqn{\texttildelow} operator, and the terms on
  the right. The response must be a survival object or Hist
  object as returned by the 'Surv' or 'Hist' function.}

  \item{data}{a data frame in which to interpret the
  variables named in the \code{formula}.}

  \item{eps}{a vector of length 3 for the convergence
  criteria (criterion for parameters, criterion for
  likelihood, criterion for second derivatives). The
  default is 'c(5,5,3)' and corresponds to criteria equals
  to \eqn{10^{-5}}, \eqn{10^{-5}} and \eqn{10^{-3}}.}

  \item{nknots}{number of knots for the splines to use to
  approximate the hazard function. Argument for the
  penalized likelihood approach.  The default is 7.}

  \item{CV}{binary variable equals to 1 when search (by
  approximated cross validation) of the smoothing parameter
  kappa and 0 otherwise. Argument for the penalized
  likelihood approach. The default is 0.}

  \item{kappa}{if CV=FALSE, smoothing parameter; if
  CV=TRUE, initial value of the smoothing parameters for
  the cross validation search. Argument for the penalized
  likelihood approach.}

  \item{conf.int}{Boolean parameter. Equals to \code{TRUE}
  to calculate pointwise confidence intervals for the
  survival or hazard curves, \code{FALSE} otherwise.
  Default is \code{TRUE}.}

  \item{maxiter}{maximum number of iterations. The default
  is 200.}

  \item{method}{type of estimation method: "Splines" for a
  penalized likelihood approach with approximation of the
  hazard function by M-splines, "Weib" for a parametric
  approach with a Weibull distribution on the hazard
  function. Default is "Weib".}

  \item{print.iter}{boolean parameter. Equals to
  \code{TRUE} to print the likelihood during the iteration
  process, \code{FALSE} otherwise. Default is \code{FALSE}.
  This option is not running on Windows.}

  \item{na.action}{how NAs are treated. The default is
  first, any na.action attribute of data, second a
  na.action setting of options, and third 'na.fail' if that
  is unset. The 'factory-fresh' default is na.omit. Another
  possible value is NULL.}
}
\value{
\item{call}{} \item{coef}{regression parameters.}
\item{loglik}{vector containing the log-likelihood without
and with covariate.} \item{modelPar}{Weibull parameters.}
\item{N}{number of subjects.} \item{NC}{number of
covariates.} \item{nevents}{number of events.}
\item{modelResponse}{model response: \code{Hist} or
\code{Surv} object.} \item{converged}{integer equal to 1
when the model converged, 2, 3 or 4 otherwise.}
\item{time}{times for which survival and hazard functions
have been evaluated for plotting.} \item{hazard}{matched
values of the hazard function.} \item{lowerHazard}{lower
confidence limits for hazard function.}
\item{upperHazard}{upper confidence limits for hazard
function.} \item{surv}{matched values of the survival
function.} \item{lowerSurv}{lower confidence limits for
survival function.} \item{upperSurv}{upper confidence
limits for survival function.} \item{RR}{vector of relative
risks.} \item{V}{variance-covariance matrix.}
\item{se}{standart errors.} \item{knots}{knots to
approximate by M-splines the hazard function.}
\item{nknots}{number of knots.} \item{CV}{a binary variable
equals to 1 when search of the smoothing parameter
\link{kappa} by approximated cross-validation, 1 otherwise.
The default is 0.} \item{niter}{number of iterations.}
\item{cv}{vector containing the convergence criteria.}
\item{na.action}{observations deleted if missing values.}
}
\description{
Fit a survival model using either a semi-parametric
approach (penalized likelihood with an approximation of the
hazard function by linear combination of M-splines) or a
parametric approach (specifying a Weibull distribution on
the hazard function). Left-truncated, right-censored, and
interval-censored data are allowed.
}
\details{
The estimated parameters are obtained using the robust
Marquardt algorithm (Marquardt, 1963) which is a
combination between a Newton-Raphson algorithm and a
steepest descent algorithm.
}
\examples{
# Weibull survival model
data(testdata)
fit.su <- shr(Hist(time=list(l,r),id)~cov,data=testdata,)

## to print
fit.su

## summary
summary(fit.su)
}
\author{
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>
}
\references{
D. Marquardt (1963). An algorithm for least-squares
estimation of nonlinear parameters.  \emph{SIAM Journal of
Applied Mathematics}, 431-441.
}
\seealso{
\code{\link{shr}}, \code{\link{print.shr}},
\code{\link{summary.shr}}, \code{\link{print.shr}},
}
\keyword{methods}

=======
\name{survival model}
\alias{shr}

\title{
Fit a survival model
}
\description{
Fit a survival model using either a semi-parametric approach (penalized likelihood with an approximation of the hazard function by linear combination of M-splines) or a parametric approach (specifying a Weibull distribution on the hazard function).
Left-truncated, right-censored, and interval-censored data are allowed.
}
\usage{
shr(formula, data, eps = c(5, 5, 3), nknots=7,
CV=FALSE , kappa = 10000, maxiter = 200, hazard = "Weib",
conf.int=TRUE, print.iter = FALSE, na.action = na.omit)
}

\arguments{
  \item{formula}{
a formula object with the response on the left of a \eqn{\texttildelow}
operator, and the terms on the right. The response must be a
survival object or Hist object as returned by the 'Surv' or 'Hist' function.
}
  \item{data}{
 a data frame in which to interpret the variables named in the \code{formula}.
}
  \item{eps}{
a vector of length 3 for the convergence criteria (criterion for parameters, criterion for likelihood, criterion for second derivatives). The default is 'c(5,5,3)' and corresponds to criteria equals to \eqn{10^{-5}}, \eqn{10^{-5}} and \eqn{10^{-3}}.
}
  \item{nknots}{
number of knots for the splines to use to approximate the hazard function. Argument for the penalized likelihood approach.
  The default is 7.
}
  \item{CV}{
binary variable equals to 1 when search (by approximated cross validation) of the smoothing parameter kappa and 0 otherwise. Argument for the penalized likelihood approach. The default is 0.
}

  \item{kappa}{
if CV=FALSE, smoothing parameter; if CV=TRUE, initial value of the smoothing parameters for the cross validation search. Argument for the penalized likelihood approach.  
}

  \item{maxiter}{
maximum number of iterations. The default is 200. 
}
  \item{hazard}{
type of estimation method: "Splines" for a penalized likelihood approach with approximation of the hazard function by M-splines, "Weib" for a parametric approach with a Weibull distribution on the hazard function. Default is "Weib".
}
  \item{conf.int}{
Boolean parameter. Equals to \code{TRUE} to calculate pointwise confidence
intervals for the survival or hazard curves, \code{FALSE} otherwise. Default is \code{TRUE}.
  }
  \item{print.iter}{
boolean parameter. Equals to \code{TRUE} to print the likelihood during the iteration process, \code{FALSE} otherwise. Default is \code{FALSE}. This option is not running on Windows.
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
  \item{call}{}
  \item{coef}{regression parameters.}
  \item{loglik}{vector containing the log-likelihood without and with covariate.} 
  \item{modelPar}{Weibull parameters.}
  \item{N}{number of subjects.}
  \item{NC}{number of covariates.}
  \item{nevents}{number of events.}
  \item{modelResponse}{model response: \code{Hist} or \code{Surv} object.} 
  \item{converged}{integer equal to 1 when the model converged, 2, 3 or 4 otherwise.}
  \item{time}{times for which survival and hazard functions have been evaluated for plotting.}
  \item{hazard}{matched values of the hazard function.}
  \item{lowerHazard}{lower confidence limits for hazard function.}  
  \item{upperHazard}{upper confidence limits for hazard function.}
  \item{surv}{matched values of the survival function.}
  \item{lowerSurv}{lower confidence limits for survival function.}
  \item{upperSurv}{upper confidence limits for survival function.} 
  \item{RR}{vector of relative risks.}
  \item{V}{variance-covariance matrix.}
  \item{se}{standart errors.}
  \item{knots}{knots to approximate by M-splines the hazard function.}
  \item{nknots}{number of knots.}
  \item{CV}{a binary variable equals to 1 when search of the smoothing parameter \link{kappa} by approximated cross-validation, 1 otherwise. The default is 0.}  
  \item{niter}{number of iterations.} 
  \item{cv}{vector containing the convergence criteria.}
  \item{na.action}{observations deleted if missing values.} 
}
\seealso{
\code{\link{shr}},
\code{\link{print.shrWeib}},
\code{\link{summary.shrWeib}},
\code{\link{print.shrSplines}},
\code{\link{summary.shrSplines}}
}
\references{
D. Marquardt (1963). An algorithm for least-squares estimation of nonlinear parameters. 
\emph{SIAM Journal of Applied Mathematics}, 431-441. 
}
\author{
R: Celia Touraine <Celia.Touraine@isped.u-bordeaux2.fr>
Fortran: Pierre Joly <Pierre.Joly@isped.u-bordeaux2.fr>
}


\examples{
\dontrun{
# Weibull survival model 
data(test)
fit.su <- shr(Hist(time=list(l,r),id)~cov,data=test) 

## to print
fit.su

## summary
summary(fit.su)
}
}

\keyword{methods}
\keyword{shrSplines}
\keyword{shrWeib}
>>>>>>> saves before pulling Celia's changes
