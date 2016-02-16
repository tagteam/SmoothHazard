#' Fit a survival model
#' 
#' Fit a survival model using either a semi-parametric approach (penalized
#' likelihood with an approximation of the hazard function by linear
#' combination of M-splines) or a parametric approach (specifying a Weibull
#' distribution on the hazard function). Left-truncated, right-censored, and
#' interval-censored data are allowed.
#' 
#' The estimated parameters are obtained using the robust Marquardt algorithm
#' (Marquardt, 1963) which is a combination between a Newton-Raphson algorithm
#' and a steepest descent algorithm.
#' 
#' @param formula a formula object with the response on the left of a
#' \eqn{\texttildelow} operator, and the terms on the right. The
#' response must be a survival object or Hist object as returned by
#' the 'Surv' or 'Hist' function.
#' @param data a data frame in which to interpret the variables named
#' in the \code{formula}.
#' @param eps a vector of length 3 for the convergence criteria
#' (criterion for parameters, criterion for likelihood, criterion for
#' second derivatives). The default is 'c(5,5,3)' and corresponds to
#' criteria equals to \eqn{10^{-5}}, \eqn{10^{-5}} and \eqn{10^{-3}}.
#' @param n.knots  Argument only active for the penalized likelihood approach \code{method="splines"}.
#' Number of knots for the splines to use to approximate
#' the hazard function. The default is 7. If \code{knots} are given as a vector this argument is ignored.
#' The algorithm needs least 5 knots and at most 20 knots.
#' @param knots Argument only active for the penalized likelihood approach \code{method="splines"}.
#' There are three ways to control the placement of the knots between the smallest and the largest
#' of all time points:
#' \itemize{
#'  \item{\code{knots="equidistant"}}{Knots are placed with same distance on the time scale.}
#'  \item{\code{knots="quantiles"}}{Knots are placed such that the number of observations is roughly the same between knots.}
#' \item{knots=list()}{List of length 3. The list elements are the actual placements
#' (timepoints) of the knots for the M-spline.}
#' }
#' The algorithm needs at least 5 knots and allows no more than 20 knots.
#' @param CV binary variable equals to 1 when search (by approximated
#' cross validation) of the smoothing parameter kappa and 0
#' otherwise. Argument for the penalized likelihood approach. The
#' default is 0.
#' @param kappa Argument only active for the penalized likelihood approach \code{method="splines"}.
#' A positive number (smoothing parameter) 
#' If CV=1 the value is used as a starting value 
#' for a cross validation search to optimize \code{kappa}.
#' @param conf.int Boolean parameter. Equals to \code{TRUE} to
#' calculate pointwise confidence intervals for the survival or hazard
#' curves, \code{FALSE} otherwise. Default is \code{TRUE}.
#' @param level confidence level for the pointwise confidence intervals 
#' of the curves. Default is 0.95.
#' @param maxiter maximum number of iterations. The default is 200.
#' @param method type of estimation method: "Splines" for a penalized
#' likelihood approach with approximation of the hazard function by
#' M-splines, "Weib" for a parametric approach with a Weibull
#' distribution on the hazard function. Default is "Weib".
#' @param print.iter boolean parameter. Equals to \code{TRUE} to print
#' the likelihood during the iteration process, \code{FALSE}
#' otherwise. Default is \code{FALSE}. This option is not running on
#' Windows.
#' @param na.action how NAs are treated. The default is first, any
#' na.action attribute of data, second a na.action setting of options,
#' and third 'na.fail' if that is unset. The 'factory-fresh' default
#' is na.omit. Another possible value is NULL.
#' @return \item{call}{} \item{coef}{regression parameters.}
#' \item{loglik}{vector containing the log-likelihood without and with
#' covariate.} \item{modelPar}{Weibull parameters.} \item{N}{number of
#' subjects.} \item{NC}{number of covariates.} \item{nevents}{number of
#' events.} \item{modelResponse}{model response: \code{Hist} or \code{Surv}
#' object.} \item{converged}{integer equal to 1 when the model converged, 2, 3
#' or 4 otherwise.} \item{time}{times for which survival and hazard functions
#' have been evaluated for plotting.} \item{hazard}{matched values of the
#' hazard function.} \item{lowerHazard}{lower confidence limits for hazard
#' function.} \item{upperHazard}{upper confidence limits for hazard function.}
#' \item{surv}{matched values of the survival function.} \item{lowerSurv}{lower
#' confidence limits for survival function.} \item{upperSurv}{upper confidence
#' limits for survival function.} \item{RR}{vector of relative risks.}
#' \item{V}{variance-covariance matrix.} \item{se}{standard errors.}
#' \item{knots}{knots of the M-splines estimate of the hazard function.}
#' \item{nknots}{number of knots.} \item{CV}{a binary variable equals to 1 when
#' search of the smoothing parameter \link{kappa} by approximated
#' cross-validation, 1 otherwise. The default is 0.} \item{niter}{number of
#' iterations.} \item{cv}{vector containing the convergence criteria.}
#' \item{na.action}{observations deleted if missing values.}
#' @author R: Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> Fortran:
#' Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' @seealso \code{\link{shr}}, \code{\link{print.shr}},
#' \code{\link{summary.shr}}, \code{\link{print.shr}},
#' @references D. Marquardt (1963). An algorithm for least-squares estimation
#' of nonlinear parameters.  \emph{SIAM Journal of Applied Mathematics},
#' 431-441.
#' @keywords methods 
##' @examples
##' 
##' # Weibull survival model
##' library(prodlim)
##' data(testdata)
##' fit.su <- shr(Hist(time=list(l,r),id)~cov,data=testdata) 
##' fit.su
##' summary(fit.su)
##'\dontrun{
##' shr.spline <- shr(Hist(time=list(l,r),id)~cov,data=testdata,method="splines",n.knots=6)
##' shr.spline
##' shr.spline.q <- shr(Hist(time=list(l,r),id)~cov,data=testdata,method="splines",n.knots=6,knots="quantiles")
##' plot(shr.spline.q)
##'
##' ## manual placement of knots
##' shr.spline.man <- shr(Hist(time=list(l,r),id)~cov,data=testdata,method="splines",knots=seq(0,7,1))
##' }
#' @export
shr <- function(formula,
                data,
                eps=c(5,5,3),
                n.knots=7,
                knots="equidistant",
                CV=FALSE,
                kappa=10000,
                conf.int=TRUE,
                level=.95,
                maxiter=200,
                method="Weib",
                print.iter=FALSE,
                na.action=na.omit){

  call <- match.call()
  ## cat("\n")
  ## cat("Be patient. The program is computing ... \n") 
	
  flush.console()
  ptm<-proc.time()
  # {{{ process formula and data
  # check if formula is a formula
  method <- tolower(method)
  if(!(method %in% c("weib","splines"))) stop("The method must be either 'Weib' or 'Splines'")
  # --------------------------------------------------------------------
  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[1]=="~")||(match("$",formula.names,nomatch=0)+match("[",formula.names,nomatch=0)>0)){
    stop("Invalid specification of formula. Perhaps forgotten right hand side?\nNote that any subsetting, ie data$var or data[,\"var\"], is invalid for this function.")
  }else{
    if (!(formula.names[2] %in% c("Surv","Hist"))) stop("formula is NOT a proper survival formula,\nwhich must have a `Surv' or `Hist' object as response.")
  }
	
  call <- match.call()
  m <- match.call()
  position <- match(c("","formula","data","na.action"),names(m),0)
  m <- m[position]
  m[[1]] <- as.name("model.frame")	
  m <- eval(m,sys.parent())
  na.action <- attr(m,"na.action")
  response <- stats::model.response(m)

  #  FIX for people who use `Surv' instead of `Hist' 
  
  if (match("Shr",class(response),nomatch=0)!=0){
    # classe "Shr"
    attr(response,"model") <- "survival"
    attr(response,"cens.type") <- "rightCensored"
    model.type <- 1
  }else{
    # classe "Hist"
    model.type <- match(attr(response,"model"),c("survival","competing.risks","multi.states"))
  }

  cens.type <- attr(response,"cens.type")
  event.history <- response

  if (cens.type!="intervalCensored"){
    event.time.order <- order(event.history[,"time"],-event.history[,"status"])
  }else{
    event.time.order <- order(event.history[,"L"],-event.history[,"status"])
  }

  # }}}
  # {{{ covariates
  # add Amadou 30/01/12
  #----------------------- factor
  X <- model.matrix(formula,data=m)
  Xnames <- NULL
  if (ncol(X) == 1){ 
    #no covariate
    X <- X-1
    NC <- 0
  }else{
    # covariate
    X <- X[, -1, drop = FALSE]
    NC <- NCOL(X)
    Xnames <- colnames(X)
  } 
  #----------------------- factor	

  # }}}
  # {{{ prepare censored event times 
  N <- NROW(X)

  isIntervalCensored <- cens.type=="intervalCensored"
	
  truncated <- nchar(attr(event.history,"entry.type"))>1  

  if (truncated){
    entrytime <- event.history[,"entry"]
  }else{
    entrytime <- rep(0,N)
  }

  if (isIntervalCensored){
    Ltime <- as.double(event.history[,"L",drop=TRUE])
    Rtime <- as.double(event.history[,"R",drop=TRUE])
  }else{
    Ltime <- as.double(event.history[,"time",drop=TRUE])
    Rtime <- rep(0,N)
  }
  id <- event.history[,"status"]
  # }}}
  # {{{ call Fortran function weib and collect results

  # do not give infinite values to fortran
  Rtime[is.infinite(Rtime)] <- Ltime[is.infinite(Rtime)]

  if (method == "weib"){	
      size1 <- NC
      size2 <- size1^2
      size_V <- size1 + 2
      ffit <- .Fortran("survWeib",
                       as.double(entrytime),
                       as.double(Ltime),
                       as.double(Rtime),
                       as.integer(id),
                       as.double(X),
                       as.integer(N),
                       as.integer(NC),
                       as.integer(truncated),
                       as.integer(isIntervalCensored),
                       as.integer(eps),
                       as.integer(maxiter),
                       loglik=as.double(rep(0,2)),
                       basepar=as.double(rep(0.1,2)),
                       regpar=as.double(rep(0.1,NC)),
                       v=as.double(rep(0,NC*NC)),
                       converged=as.integer(rep(0,2)),
                       cv=as.double(rep(0,3)),
                       niter=as.integer(0),
                       t=as.double(rep(0,100)),
                       S=as.double(rep(0,100)),
                       S_l=as.double(rep(0,100)),
                       S_u=as.double(rep(0,100)),
                       h=as.double(rep(0,100)),
                       h_l=as.double(rep(0,100)),
                       h_u=as.double(rep(0,100)),
                       as.integer(conf.int),
                       as.double(level),
                       as.integer(print.iter),
                       V_tot=as.double(matrix(0,nrow=size_V,ncol=size_V)),
                       PACKAGE="SmoothHazard")
  }else{
       #  	cat("------ Program Splines ------ \n")
       if (length(entrytime)>0){
           alltimes <- sort(unique(c(Ltime, Rtime,entrytime)))
           amax <- max(alltimes)
           amin <- min(alltimes)
       }
       else{
           alltimes <- sort(unique(c(Ltime, Rtime)))
           amax <- max(alltimes)
           amin <- min(alltimes)
       }
       if (is.character(knots)){
           if (length(n.knots)!=1) stop("Argument n.knots has to be a single positive integer")
           if((!is.numeric(n.knots) && !is.integer(n.knots)) || (n.knots < 5) || (n.knots >20))
               stop("n.knots has to be an integer between 5 and 20. See help(shr).")
           if (knots=="quantiles"){
               knots <- quantile(alltimes,seq(0,1,1/(n.knots-1)))
           } else {
                 if (knots!="equidistant")
                     warning("Unknown specification of knots. Fall back to equidistant.")
                 knots <- seq(amin,amax,(amax-amin)/(n.knots-1))
             }
       } else{## user specified knots
             if (!is.numeric(knots))
                 stop("Incorrect form of argument knots. See help(shr).")
             if (knots[1]< amin - 0.05*amin) stop(paste("Smallest knot should not be smaller than the time point:",amin))
             if (knots[length(knots)]> amax + 0.05*amax) stop(paste("Largest knot should not be larger than the time point:",amax))
             knots <- sort(knots)
         }
       nknots <- length(knots)
       if (nknots<5) {
           stop("Need at least 5 knots.")
       }
       if (nknots>20){
           stop("Cannot handle more than 20 knots.")
       }
       ## make fake knots needed for M-splines
       knots <- c(rep(knots[1],3),knots,rep(knots[length(knots)],3))
       size1 <- NC
       size2 <- size1^2
       size_V <- size1 + nknots + 2
       ffit <- .Fortran("survPl",
                        as.double(entrytime),
                        as.double(Ltime),
                        as.double(Rtime),
                        as.integer(id),
                        as.double(X),
                        as.integer(N),
                        as.integer(NC),
                        as.integer(truncated),
                        as.integer(isIntervalCensored),
                        as.integer(eps),
                        as.integer(maxiter),
                        loglik=as.double(rep(0,2)),
                        regpar=as.double(rep(0.1,NC)),
                        v=as.double(rep(0,NC*NC)),
                        converged=as.integer(rep(0,2)),
                        cv=as.double(rep(0,3)),
                        niter=as.integer(0),
                        t=as.double(rep(0,99)),
                        S=as.double(rep(0,99)),
                        S_l=as.double(rep(0,99)),
                        S_u=as.double(rep(0,99)),
                        h=as.double(rep(0,99)),
                        h_l=as.double(rep(0,99)),
                        h_u=as.double(rep(0,99)),
                        as.integer(nknots),
                        as.double(knots),
                        as.integer(CV),
                        as.double(kappa),
                        kappa=as.double(0),
                        as.integer(conf.int),
                        as.double(level),
                        CVcrit=as.double(0),
                        mdf=as.double(0),
                        ti=as.double(rep(0,(nknots+6))),
                        theta=as.double(rep(0,(nknots+2))),
                        as.integer(print.iter),
                        V_tot=as.double(matrix(0,nrow=size_V,ncol=size_V)),
                        PACKAGE="SmoothHazard")



  }
  ## Fortran delivers
  ## Variable name 	Explanation 	Dimension 	Storage mode 	Remark
  ## loglik 	log-likelihood without and with covariate 	length 2 	double 	
  ## basepar 	Weibull parameters 	length 2 	double 	
  ## regpar 	Regression coefficients 	length P 	double 	
  ## v 	covariance matrix 	length P*P 	double 	
  ## converged 	0=converged,1=invert fails,2=no 	length 1 	integer 	
  ## t 	time to plot S(t) and h(t) 	length 100 	double 	
  ## S 	survival function 	length 100 	double 	
  ## S_l 	lower confidence limit for S 	length 100 	double 	
  ## S_u 	Upper confidence limit for S 	length 100 	double 	
  ## h 	hazard function 	length 100 	double 	
  ## h_l 	lower confidence limit for h function 	length 100 	double 	
  ## h_u 	upper confidence limit for h 	length 100 	double 	

  if (any(ffit$converged == 4)){
      warning("Problem in the loglikelihood computation. The program stopped abnormally. Please check your dataset. \n")    
  }
  if (any(ffit$converged == 2)){
      if (CV==0) 
          warning("Model did not converge. You could try to increase the 'maxit' parameter and set 'CV=1'.")
      else
          warning("Model did not converge. You could try to increase the 'maxit' parameter and modify the start values for 'kappa'.")
  }
  if (any(ffit$converged == 3)){
      warning("The Fisher information matrix is not positive definite.")
  }
  if(method=="weib") weibullParameter <- ffit$basepar
  # }}}
  # {{{ output 
	
  fit <- NULL
  fit$call <- call
	

  fit$loglik <- ffit$loglik
  if(method=="weib"){
    fit$modelPar <- weibullParameter
  }
  fit$N <- N
  fit$NC <- NC
  fit$modelResponse <- event.history
  fit$converged <- ffit$converged
  fit$time <- ffit$t
  fit$hazard <- ffit$h
  fit$lowerHazard <- ffit$h_l
  fit$upperHazard <- ffit$h_u
  fit$surv <- ffit$S
  fit$lowerSurv <- ffit$S_l
  fit$upperSurv <- ffit$S_u
	
  if(NC>0){
    betaCoef <- ffit$regpar
    names(betaCoef) <- Xnames
    fit$coef <- betaCoef
    fit$HR <- exp(betaCoef)
    V <- matrix(ffit$v,nrow=NC,ncol=NC,byrow=T)
    colnames(V) <- Xnames
    rownames(V) <- Xnames
    fit$se <- sqrt(diag(V))
    fit$V_cov <- V
  }
  V <- matrix(ffit$V_tot,nrow=size_V,ncol=size_V,byrow=T)
  if(method=="weib"){
    colnames(V) <- c("sqrt(a)","sqrt(b)",Xnames)
    rownames(V) <- c("sqrt(a)","sqrt(b)",Xnames)
  }else{
    theta_names <- cbind(rep("theta",(nknots+2)),(1:(nknots+2)))
    theta_names <- as.vector(apply(theta_names,1,paste,collapse=" "))
    colnames(V) <- c(theta_names,Xnames)	
    rownames(V) <- c(theta_names,Xnames)	
  }
  fit$V <- V
  fit$niter <- ffit$niter
  fit$cv <- ffit$cv

  if(method=="splines"){
    fit$nknots <- nknots
    fit$knots <- ffit$ti
    fit$theta <-  ffit$theta
    fit$CV <- CV
    fit$igraph <- ffit$igraph
    if(CV){
      fit$kappa <- ffit$kappa
      fit$CVcrit <- ffit$CVcrit
      fit$DoF <- ffit$mdf
    }else{
      fit$kappa <- kappa
    }

}
  fit$na.action <- na.action
  # }}}
  if (method=="weib") fit$method <- "Weib" else fit$method <- "Splines"
  class(fit) <- "shr"
  fit$runtime <- proc.time()-ptm
  fit
}
