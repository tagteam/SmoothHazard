#' Fit an illness-death model
#' 
#' Fit an illness-death model using either a semi-parametric approach
#' (penalized likelihood with an approximation of the transition intensity
#' functions by linear combination of M-splines) or a parametric approach
#' (specifying Weibull distributions on the transition intensities).
#' Left-truncated, right-censored, and interval-censored data are allowed.
#' State 0 corresponds to the initial state, state 1 to the transient one,
#' state 2 to the absorbant one. The allowed transitions are: 0 --> 1, 0 --> 2
#' and 1 --> 2.
#' 
#' The estimated parameters are obtained using the robust Marquardt algorithm
#' (Marquardt, 1963) which is a combination between a Newton-Raphson algorithm
#' and a steepest descent algorithm.
#' 
#' @param formula01 A formula specifying a regression model for the
#' \code{0 --> 1} transition from the initial state to the transient
#' state of the illness-death model.  The right hand side of the
#' formula specifies the covariate terms, and the left hand side must
#' be an event history object as returned by the function \code{Hist}.
#' @param formula02 A formula specifying a regression model for the
#' \code{0 --> 2} transition from the initial state to the absorbing
#' state. The left hand side must be equal to the left hand side of
#' \code{formula01}. If missing it is set to \code{formula01}.
#' @param formula12 A formula specifying a regression model for the
#' \code{1 --> 2} transition from the transient state to the absorbing
#' state.  operator is not required. If missing it is set to
#' \code{formula01}.
#' @param data A data frame in which to interpret the variables of
#' \code{formula01}, \code{formula02} and \code{formula12}.
#' @param maxiter Maximum number of iterations. The default is 200.
#' @param eps A vector of 3 integers >0 used to define the power of
#' three convergence criteria: 1. for the regression parameters,
#' 2. for the likelihood, 3. for the second derivatives. The default
#' is \code{c(5,5,3)} which is translated into convergence if the
#' respective values change less then \eqn{10^{-5}} (for regression
#' parameters and likelihood) and \eqn{10^{-3}} for the second
#' derivatives between two iterations.
#' @param n.knots For \code{method="Splines"} only, a vector of length
#' 3 specifing the number of knots, one for each transition, for the
#' M-splines estimate of the baseline intensities in the order \code{0
#' --> 1}, \code{0 --> 2}, \code{1 --> 2}. The default is c(7,7,7). When \code{knots}
#' are specified as a list this argument is ignored.
#' The algorithm needs least 5 knots and at most 20 knots.
#' @param knots Argument only active for the penalized likelihood approach \code{method="splines"}.
#' There are three ways to control the placement of the knots between the smallest and the largest
#' of all time points:
#' \itemize{
#'  \item{\code{knots="equidistant"}}{Knots are placed with same distance on the time scale.}
#'  \item{\code{knots="quantiles"}}{Knots are placed such that the number of observations is roughly the same between knots.}
#' \item{knots=list()}{List of 1 or 2 or three vectors. The list elements are the actual placements
#' (timepoints) of the knots for the M-spline. The list may contain 
#' one vector of placements for each transition in the order \code{0 --> 1}, \code{0 --> 2}, \code{1 --> 2}.
#' If only vector is specifified the knots are used for all transitions. If only 2 vectors are specifified, the
#' knots for the \code{0 --> 1} transition are also used for the \code{1 --> 2} transition.}
#' }
#' The algorithm needs at least 5 knots and allows no more than 20 knots.
#' @param CV Binary variable equals to 1 when search (by approximated
#' cross validation) of the smoothing parameters \code{kappa} and 0
#' otherwise. Argument for the penalized likelihood approach. The
#' default is 0.
#' @param kappa Argument only active for the penalized likelihood approach \code{method="splines"}.
#' A vector with 3 positive values (smoothing parameters), one for each transition, in the order 
#' 0 --> 1, 0 --> 2 and 1 --> 2..
#' If CV=1 these are used as starting values for a cross validation search to optimize kappa.
#' @param method type of estimation method: "Splines" for a
#' penalized likelihood approach with approximation of the transition
#' intensities by M-splines, "Weib" for a parametric approach with a
#' Weibull distribution on the transition intensities. Default is
#' "Weib".
#' @param conf.int Logical. If \code{TRUE} 
#' calculate pointwise confidence intervals for the transition
#' intensities curves, \code{FALSE} otherwise. Default is \code{TRUE}.
#' @param level confidence level for the pointwise confidence intervals 
#' of the transition intensities curves. Default is 0.95.
#' @param print.iter boolean parameter. Equals to \code{TRUE} to print
#' the likelihood during the iteration process, \code{FALSE}
#' otherwise. Default is \code{FALSE}. This option is not running on
#' Windows.
#' @param subset expression indicating the subset of the rows of data
#' to be used in the fit. All observations are included by default.
#' @param na.action how NAs are treated. The default is first, any
#' na.action attribute of data, second a na.action setting of options,
#' and third 'na.fail' if that is unset. The 'factory-fresh' default
#' is na.omit. Another possible value is NULL.
#' @return
#' 
#' \item{call}{the call that produced the result.} \item{coef}{regression
#' parameters.} \item{loglik}{vector containing the log-likelihood without and
#' with covariate.} \item{cv}{vector containing the convergence criteria.}
#' \item{niter}{number of iterations.} \item{converged}{integer equal to 1 when
#' the model converged, 2, 3 or 4 otherwise.} \item{modelPar}{Weibull
#' parameters.} \item{N}{number of subjects.} \item{events1}{number of events 0
#' --> 1.} \item{events2}{number of events 0 --> 2 or 0 --> 1 --> 2.}
#' \item{NC}{vector containing the number of covariates on transitions 0 --> 1,
#' 0 --> 2, 1 --> 2.} \item{responseTrans}{model response for the 0 --> 1
#' transition. \code{Hist} or \code{Surv} object.} \item{responseAbs}{model
#' response for the 0 --> 2 transition. \code{Hist} or \code{Surv} object.}
#' \item{time}{times for which transition intensities have been evaluated for
#' plotting. Vector in the Weibull approach. Matrix in the penalized likelihhod
#' approach for which the colums corresponds to the transitions 0 --> 1, 1 -->
#' 2, 0 --> 2.} \item{intensity01}{matched values of the intensities for
#' transition 0 --> 1.} \item{lowerIntensity01}{lower confidence intervals for
#' the values of the intensities for transition 0 --> 1.}
#' \item{upperIntensity01}{upper confidence intervals for the values of the
#' intensities for transition 0 --> 1.} \item{intensity02}{matched values of
#' the intensities for transition 0 --> 2.} \item{lowerIntensity02}{lower
#' confidence intervals for the values of the intensities for transition 0 -->
#' 2.} \item{upperIntensity02}{upper confidence intervals for the values of the
#' intensities for transition 0 --> 2.} \item{intensity12}{matched values of
#' the intensities for transition 1 --> 2.} \item{lowerIntensity12}{lower
#' confidence intervals for the values of the intensities for transition 1 -->
#' 2.} \item{upperIntensity12}{upper confidence intervals for the values of the
#' intensities for transition 1 --> 2.} \item{RR}{vector of relative risks.}
#' \item{V}{variance-covariance matrix.} \item{se}{standart errors of the
#' regression parameters.} \item{Xnames01}{names of covariates on 0 --> 1.}
#' \item{Xnames02}{names of covariates on 0 --> 2.} \item{Xnames12}{names of
#' covariates on 1 --> 2.} \item{knots01}{knots to approximate by M-splines the
#' intensity of the 0 --> 1 transition.} \item{knots02}{knots to approximate by
#' M-splines the intensity of the 0 --> 2 transition.} \item{knots12}{knots to
#' approximate by M-splines the intensity of the 1 --> 2 transition.}
#' \item{nknots01}{number of knots on transition 0 --> 1.}
#' \item{nknots02}{number of knots on transition 0 --> 2.}
#' \item{nknots12}{number of knots on transition 1 --> 2.}
#' \item{theta01}{square root of splines coefficients for transition 0 --> 1.}
#' \item{theta02}{square root of splines coefficients for transition 0 --> 2.}
#' \item{theta12}{square root of splines coefficients for transition 1 --> 2.}
#' \item{CV}{a binary variable equals to 1 when search of the smoothing
#' parameters \link{kappa} by approximated cross-validation, 1 otherwise. The
#' default is 0.} \item{kappa}{vector containing the smoothing parameters for
#' transition 0 --> 1, 0 --> 2, 1 --> 2 used to estimate the model by the
#' penalized likelihood approach.} \item{CVcrit}{cross validation criteria.}
#' \item{DoF}{degrees of freedom of the model.} \item{na.action}{observations
#' deleted if missing values.}
#' @author R: Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> Fortran:
#' Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' @seealso \code{\link{print.idm}}
#' \code{\link{summary.idm}}
#' \code{\link{predict.idm}}
#' @references D. Marquardt (1963). An algorithm for least-squares estimation
#' of nonlinear parameters.  \emph{SIAM Journal of Applied Mathematics},
#' 431-441.
#' @keywords illness-death
#' 
##' @examples
##' library(lava)
##' library(prodlim)
##' set.seed(17)
##' d <- simulateIDM(100)
##' # right censored data
##' fitRC <- idm(formula01=Hist(time=observed.illtime,event=seen.ill)~X1+X2,
##'              formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
##'              formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d,
##'              conf.int=FALSE)
##' fitRC
##'
##' \dontrun{
##' fitRC.splines <- idm(formula01=Hist(time=observed.illtime,event=seen.ill)~X1+X2,
##'              formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
##'              formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d,
##'              conf.int=FALSE,method="splines")
##' }
##' # interval censored data
##' fitIC <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
##'              formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
##'              formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d,
##'              conf.int=FALSE)
##' fitIC
##' 
##' \dontrun{
##' 
##'     data(Paq1000)
##' 
##'     # Illness-death model with certif on the 3 transitions
##'     # Weibull parametrization and likelihood maximization
##' 
##'     fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
##'                     formula01=Hist(time=list(l,r),event=dementia)~certif,
##'                     data=Paq1000)
##' 
##'     # Illness-death model with certif on transitions 01 and 02
##'     # Splines parametrization and penalized likelihood maximization
##'     fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=e)~certif,
##'                         formula01=Hist(time=list(l,r),event=dementia)~certif,
##'                         formula12=~1,
##'                         method="Splines",
##'                         data=Paq1000)
##'     fit.weib
##'     summary(fit.splines)
##' }
##' 
#' @importFrom prodlim Hist
#' @useDynLib SmoothHazard
#' @export
idm <- function(formula01,
                formula02,
                formula12,
                data,
                maxiter=200,
                eps=c(5,5,3),
                n.knots=c(7,7,7),
                knots="equidistant",
                CV=FALSE,
                kappa=c(1000000,500000,20000),
                method="Weib",
                conf.int=TRUE,
                level=.95,
                print.iter=FALSE,
                subset=NULL,
                na.action = na.fail){
    # {{{ check formula
    call <- match.call()
    ptm <- proc.time()
    if(missing(formula01))stop("Argument formula01 is missing.")
    if(missing(formula02))stop("Argument formula02 is missing.")	
    if(class(formula01)!="formula")stop("The argument formula01 must be a class 'formula'.")	
    if(class(formula02)!="formula")stop("The argument formula02 must be a class 'formula'.")		
    ## if(missing(formula02)) formula02 <- formula01
    if(missing(formula12)) formula12 <- formula02
    # }}}
    # {{{ evaluate formula in data 
    if(missing(data)) stop("Need a data frame.")
    if(class(data)!="data.frame")stop("Argument 'data' must be a data.frame")
    m <- match.call()
    m01 <- m02 <- m12 <- m[match(c("","data","subset","na.action"),names(m),nomatch=0)]
    m01$formula <- formula01
    m02$formula <- formula02
    m12$formula <- formula12
    m01[[1]] <- m02[[1]] <- m12[[1]] <- as.name("model.frame")
    m01 <- eval(m01,parent.frame())
    m02 <- eval(m02,parent.frame())
    m12 <- eval(m12,parent.frame())
    # }}}
    # {{{ check parameters
    if (any(kappa<=0)) stop("Parameter kappa has to be positive.")
    # }}}
    # {{{ extract response
    responseTrans <- stats::model.response(m01)
    responseAbs <- stats::model.response(m02)
    # }}}
    # {{{ extract covariates
    ## formula01
    x01 <- model.matrix(formula01,data=m01)[, -1, drop = FALSE]
    NC01 <- NCOL(x01)
    if (NC01>0)
        Xnames01 <- colnames(x01)
    else 
        Xnames01 <- NULL
    ## formula02
    x02 <- model.matrix(formula02,data=m02)[, -1, drop = FALSE]
    NC02 <- NCOL(x02)
    if (NC01>0)
        Xnames02 <- colnames(x02)
    else
        Xnames02 <- NULL
    ## formula12
    x12 <- model.matrix(formula12,data=m12)[, -1, drop = FALSE]
    NC12 <- NCOL(x12)
    if (NC12>0)
        Xnames12 <- colnames(x12)
    else
        Xnames12 <- NULL
    # }}}
    # {{{ prepare censored event times
    isIntervalCensored <- attr(responseTrans,"cens.type")=="intervalCensored"
    truncated <- nchar(attr(responseAbs,"entry.type"))>1
    abstime <- as.double(responseAbs[,"time"])
    ## It may happen that the illness time is observed exactly, in which case
    ## the status is 1, thus we need two criteria to declare illness status:
    ## 1. exact observations with illness status ==1
    ## 2. interval censored with any illness status. FIXME: check the corresponding likelihood
    idm <- responseTrans[,"status"]==(as.integer(isIntervalCensored)+1)
    if (isIntervalCensored)
        idm[(responseTrans[,"status"]==1 & (responseTrans[,"L"]==responseTrans[,"R"]))] <- 1
    ## exit status
    idd <- responseAbs[,"status"]==1
    N <- length(abstime)
    if (truncated==0){
        entrytime <- as.double(NULL)
    }else{
         entrytime <- as.double(responseAbs[,"entry"])
     }
    if (isIntervalCensored){
        Ltime <- as.double(responseTrans[,"L",drop=TRUE])
        Rtime <- as.double(responseTrans[,"R",drop=TRUE])
        ## if (any(Rtime<abstime & idm ==0))
        ## warning(paste("For ",
        ## sum(Rtime<abstime & idm ==0),
        ## " cases where the ill status is not observed\n and the last inspection time (R) is smaller than the right censored time (T)\n the time R is set to T."))
    }else{# exactly observed transition times
         Ltime <- as.double(responseTrans[,"time",drop=TRUE])
         Rtime <- as.double(responseTrans[,"time",drop=TRUE])
         Ltime[idm==0] <- abstime[idm==0]
         Rtime[idm==0] <- abstime[idm==0]
     }
    ## print(head(cbind(Ltime,Rtime)))
    method <- tolower(method)
    if(!(method %in% c("weib","splines"))) stop("The method argument must be 'Weib' or 'Splines'")
    # }}}
    # {{{ check data for integrity
    if (attr(responseAbs,"cens.type")=="intervalCensored") stop("No method available when the transtion to the absorbing state is interval censored.")
    if (isIntervalCensored && any(Rtime<Ltime)) stop("Misspecified transitition times:\nSome left interval limits are greater than the corresponding right limits.")
    # }}}
    # {{{ call Fortran function weib and collect results
    ## ===R provides===
    ## Variable name| Explanation|Dimension|Storage mode|Remark
    ## entrytime| truncation time|length N|double|if is_truncated = 0 then length 0
    ## l| right censored  time or left border of interval censored obs for illness|length N|double|
    ## r| right border of interval censored obs for illness|length N|double|without interval censoring r=l 
    ## d| right censored  time or obs time for death|length N|double|
    ## idm| illness status (1= event, 0=no event)|length N|integer|
    ## idd| death status (1= event, 0=no event)|length N|integer|
    ## x01| covariate matrix in vector form for 0---> 1: c(X_11,X_12,...,X_1P,X_21,X22,...,X_P*N)|length N*P|double|dummy variables
    ## x02| covariate matrix in vector form for 0---> 2: c(X_11,X_12,...,X_1P,X_21,X22,...,X_P*N)|length N*P|double|
    ## x12| covariate matrix in vector form for 1---> 2: c(X_11,X_12,...,X_1P,X_21,X22,...,X_P*N)|length N*P|double|
    ## cluster| NOT YET |--|--|
    ## strata| NOT YET |--|--|
    ## N| number of subjects|length 1|integer|
    ## P01| number of covariates for 0---> 1|length 1|integer|
    ## P02| number of covariates for 0---> 2|length 1|integer|
    ## P12| number of covariates for 1---> 2|length 1|integer|
    ## is_truncated|0=no, 1=yes|length 1|integer|
    ## eps|convergence criteria: 1:likelihood,2:parameter est,3:gradient parameter est |length 3|integer|example eps=c(7,4,5) then use 10^-7,10^-4,10^-5. Defaults to c(5,5,3)
    ## maxiter| maximum number of iteration | length 1 | integer | > 0 default to 200
    if (method == "weib"){
        #	cat("------ Program Weibull ------ \n")
        size1 <- NC01 + NC02 + NC12
        size2 <- size1^2
        size_V <- size1 + 6
        ffit <- .Fortran("idmWeib",
                         ## input
                         as.double(entrytime),               #
                         as.double(Ltime),                   #l=
                         as.double(Rtime),                   #r=
                         as.double(abstime),                 #d=
                         as.integer(idm),                    #idm=}
                         as.integer(idd),                    #idd=
                         as.double(x01),                     #x01=
                         as.double(x02),                     #x02=
                         as.double(x12),                     #x12=
                         as.integer(N),                      #N=
                         as.integer(NC01),                   #P01= 
                         as.integer(NC02),                   #P02= 
                         as.integer(NC12),                   #P12= 
                         as.integer(truncated),              #truncated=
                         ## interval=as.integer(isIntervalCensored),
                         as.integer(eps),   #eps=
                         as.integer(maxiter),
                         ## output
                         loglik=as.double(rep(0,2)),
                         basepar=as.double(rep(0,6)),
                         regpar=as.double(rep(0,size1)),
                         v=as.double(rep(0,size2)),
                         converged=as.integer(rep(0,2)),
                         cv=as.double(rep(0,3)),
                         niter=as.integer(0),
                         t=as.double(rep(0,99)),
                         a01=as.double(rep(0,99)),
                         a01_l=as.double(rep(0,99)),
                         a01_u=as.double(rep(0,99)),
                         a02=as.double(rep(0,99)),
                         a02_l=as.double(rep(0,99)),
                         a02_u=as.double(rep(0,99)),
                         a12=as.double(rep(0,99)),
                         a12_l=as.double(rep(0,99)),
                         a12_u=as.double(rep(0,99)),
                         as.integer(conf.int),
                         as.double(level),
                         as.integer(print.iter),
                         V_tot=as.double(matrix(0,nrow=size_V,ncol=size_V)),
                         PACKAGE="SmoothHazard")
    }else{
         #  	cat("------ Program Splines ------ \n")
         ## check knots
         if (length(entrytime)>0){
             alltimes <- sort(unique(c(Ltime, Rtime,entrytime,abstime)))
             amax <- max(alltimes)
             amin <- min(alltimes)
         }
         else{
             alltimes <- sort(unique(c(Ltime, Rtime,abstime)))
             amax <- max(alltimes)
             amin <- min(alltimes)
         }
         if (is.character(knots)){
             if ((length(n.knots)>3) || (length(n.knots)<1)) stop("Argument n.knots has to be a vector of at least one positive integer and at most 3 positive integers.")
             if (length(n.knots)==1) n.knots <- c(n.knots,n.knots,n.knots)
             if (length(n.knots)==2) n.knots <- c(n.knots,n.knots[1])
             nknots01 <- n.knots[1]
             nknots02 <- n.knots[2]
             nknots12 <- n.knots[3]
             if((!is.numeric(n.knots) && !is.integer(n.knots)) || (any(n.knots < 5)) || (any(n.knots >20)))
                 stop("Each element of n.knots has to be an integer between 5 and 20. See help(idm).")
             if (knots=="quantiles"){
                 approx.illtimes <- (Rtime[idm==1] + Ltime[idm==1])/2
                 knots01 <- quantile(approx.illtimes,seq(0,1,1/(nknots01-1)))
                 knots02 <- quantile(abstime,seq(0,1,1/(nknots02-1)))
                 knots12 <- quantile(abstime,seq(0,1,1/(nknots12-1)))
             }
             if (knots!="equidistant")
                 warning("Unknown specification of knots. Fall back to equidistant.")
             knots01 <- seq(amin,amax,(amax-amin)/(nknots01-1))
             knots02 <- seq(amin,amax,(amax-amin)/(nknots02-1))
             knots12 <- seq(amin,amax,(amax-amin)/(nknots12-1))
         } else{## user specified knots
               if (!is.list(knots) || length(knots)==1)
                   knots <- list(knots,knots,knots)
               if (length(knots)==2) ## re-use knots from 0->1 for 1->2
                   knots <- c(knots,knots[1])
               if (!all(sapply(knots,is.numeric)))
                   stop("Incorrect form of argument knots. See help(idm).")
               knots01 <- sort(knots[[1]])
               knots02 <- sort(knots[[2]])
               knots12 <- sort(knots[[3]])
               if (knots01[1]< amin - 0.05*amin) stop(paste("Transition 0->1: Smallest knot should not be smaller than the time point:",amin))
               if (knots01[length(knots01)]> amax + 0.05*amax) stop(paste("Transition 0->1: Largest knot should not be larger than the time point:",amax))
               if (knots02[1]< amin - 0.05*amin) stop(paste("Transition 0->2: Smallest knot should not be smaller than the time point:",amin))
               if (knots02[length(knots02)]> amax + 0.05*amax) stop(paste("Transition 0->2: Largest knot should not be larger than the time point:",amax))
               if (knots12[1]< amin - 0.05*amin) stop(paste("Transition 1->2: Smallest knot should not be smaller than the time point:",amin))
               if (knots12[length(knots12)]> amax + 0.05*amax) stop(paste("Transition 1->2: Largest knot should not be larger than the time point:",amax))
               ## FIXME: check if knots within amin, amax
               ## if (knots01[[1]] < amin) stop("Smallest knot ")
               nknots01 <- length(knots01)
               nknots02 <- length(knots02)
               nknots12 <- length(knots12)
           }
         ## double check to avoid crash
         if (any(c(nknots01,nknots02,nknots12)<5)) {
             stop("Need at least 5 knots.")
         }
         if (any(c(nknots01,nknots02,nknots12)>20)){
             stop("Cannot handle more than 20 knots.")
         }
         ## make sure min and max times are knots
         if (min(knots01)>amin) knots01 <- c(amin,knots01)
         ## warning("The first knot for the 0->1 transition has to be before or at the smallest entrytime")
         if (min(knots02)>amin) knots02 <- c(amin,knots02)
         ## stop("The first knot for the 0->2 transition has to be before or at the smallest entrytime")
         if (min(knots12)>min(Ltime)) knots02 <- c(min(Ltime,knots02))
         ## stop("The first knot for the 1->2 transition has to be before or at the first (in time) observation in the ill-state.")
         ## make fake knots needed for M-splines
         knots01 <- c(rep(knots01[1],3),knots01,rep(knots01[length(knots01)],3))
         knots02 <- c(rep(knots02[1],3),knots02,rep(knots02[length(knots02)],3))
         knots12 <- c(rep(knots12[1],3),knots12,rep(knots12[length(knots12)],3))
         size1 <- NC01 + NC02 + NC12
         size_V <- size1 + nknots01+nknots02+nknots12 + 6
         size2 <- size1**2
         ffit <- .Fortran("idmPl",
                          ## input
                          as.double(entrytime),               #entrytime=
                          as.double(Ltime),                   #l=
                          as.double(Rtime),                   #r=
                          as.double(abstime),                 #d=
                          as.integer(idm),                    #idm=
                          as.integer(idd),                    #idd=
                          as.double(x01),                     #x01=
                          as.double(x02),                     #x02=
                          as.double(x12),                     #x12=
                          as.integer(N),                      #N
                          as.integer(NC01),                   #P01= 
                          as.integer(NC02),                   #P02= 
                          as.integer(NC12),                   #P12= 
                          as.integer(truncated),              #truncated=
                          ## interval=as.integer(isIntervalCensored),
                          as.integer(eps),   #eps=
                          as.integer(maxiter),
                          ## output
                          loglik=as.double(rep(0,2)),
                          regpar=as.double(rep(0,size1)),
                          v=as.double(rep(0,size2)),
                          converged=as.integer(rep(0,2)),
                          cv=as.double(rep(0,3)),
                          niter=as.integer(0),
                          t=as.double(matrix(0,nrow=99,ncol=3)),
                          a01=as.double(rep(0,99)),
                          a01_l=as.double(rep(0,99)),
                          a01_u=as.double(rep(0,99)),
                          a02=as.double(rep(0,99)),
                          a02_l=as.double(rep(0,99)),
                          a02_u=as.double(rep(0,99)),
                          a12=as.double(rep(0,99)),
                          a12_l=as.double(rep(0,99)),
                          a12_u=as.double(rep(0,99)),	
                          as.integer(nknots01),
                          as.double(knots01),
                          as.integer(nknots02),
                          as.double(knots02),
                          as.integer(nknots12),
                          as.double(knots12),
                          as.integer(CV),
                          as.double(kappa),
                          kappa=as.double(rep(0,3)),
                          as.integer(conf.int),
                          as.double(level),
                          CVcrit=as.double(0),
                          mdf=as.double(0),
                          theta01=as.double(rep(0,(nknots01+2))),
                          theta02=as.double(rep(0,(nknots02+2))),
                          theta12=as.double(rep(0,(nknots12+2))),
                          as.integer(print.iter),
                          V_tot=as.double(matrix(0,nrow=size_V,ncol=size_V)),
                          PACKAGE="SmoothHazard")
     }
    # ffit$converged[[1]]  without covariates
    # ffit$converged[[2]]  with covariates
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
    fit <- NULL
    if(method=="weib"){
        weibullParameter <- ffit$basepar
    }
    NC <- c(NC01,NC02,NC12)
    fit$call <- call
    fit$terms <- list("Formula01"=terms(formula01),
                      "Formula02"=terms(formula02),
                      "Formula12"=terms(formula12))
    fit$loglik <- ffit$loglik
    fit$cv <- ffit$cv
    fit$niter <- ffit$niter
    fit$converged <- ffit$converged
    if(method=="weib"){
        fit$modelPar <- weibullParameter
    }
    fit$N <- N
    fit$events1 <- sum(idm)
    fit$events2 <- sum(idd)
    fit$NC <- NC
    fit$responseAbs <- responseAbs
    fit$responseTrans <- responseTrans
    if(method=="splines"){
        fit$time <- matrix(ffit$t,ncol=3) 
    }else{
         fit$time <- ffit$t 
     }
    fit$intensity01 <- ffit$a01
    fit$lowerIntensity01 <- ffit$a01_l
    fit$upperIntensity01 <- ffit$a01_u
    #
    fit$intensity02 <- ffit$a02
    fit$lowerIntensity02 <- ffit$a02_l
    fit$upperIntensity02 <- ffit$a02_u
    #
    fit$intensity12 <- ffit$a12
    fit$lowerIntensity12 <- ffit$a12_l
    fit$upperIntensity12 <- ffit$a12_u
    if (sum(NC)>0){                       # if at least one covariate
        betaCoef <- ffit$regpar
        names(betaCoef) <- c(Xnames01,Xnames02,Xnames12)
        fit$coef <- betaCoef
        fit$HR <- exp(betaCoef)
        V <- matrix(ffit$v,nrow=size1,ncol=size1,byrow=T) 
        colnames(V) <- c(Xnames01,Xnames02,Xnames12)
        rownames(V) <- c(Xnames01,Xnames02,Xnames12)
        fit$V_cov <- V
        fit$se <- sqrt(diag(fit$V))
    }  	
    V <- matrix(ffit$V_tot,nrow=size_V,ncol=size_V,byrow=T)
    if(method=="weib"){
        colnames(V) <- c("sqrt(a01)","sqrt(b01)","sqrt(a02)","sqrt(b02)","sqrt(a12)","sqrt(b12)",c(Xnames01,Xnames02,Xnames12))
        rownames(V) <- c("sqrt(a01)","sqrt(b01)","sqrt(a02)","sqrt(b02)","sqrt(a12)","sqrt(b12)",c(Xnames01,Xnames02,Xnames12))
    }else{
         theta_names <- cbind(c(rep("theta01",(nknots01+2)),rep("theta02",(nknots02+2)),rep("theta12",(nknots12+2))),c((1:(nknots01+2)),(1:(nknots02+2)),(1:(nknots12+2))))
         theta_names <- as.vector(apply(theta_names,1,paste,collapse=" "))
         colnames(V) <- c(theta_names,c(Xnames01,Xnames02,Xnames12))	
         rownames(V) <- c(theta_names,c(Xnames01,Xnames02,Xnames12))	
     }
    fit$V <- V
    if(NC01>0) fit$Xnames01 <- Xnames01
    if(NC02>0) fit$Xnames02 <- Xnames02
    if(NC12>0) fit$Xnames12 <- Xnames12
    if(method=="splines"){
        fit$knots01 <- knots01
        fit$knots02 <- knots02
        fit$knots12 <- knots12
        fit$theta01 <- ffit$theta01     
        fit$theta02 <- ffit$theta02    
        fit$theta12 <- ffit$theta12 
        fit$CV <- CV
        fit$nknots01 <- nknots01
        fit$nknots02 <- nknots02
        fit$nknots12 <- nknots12
        fit$igraph <- ffit$igraph
        fit$CVcrit <- ffit$CVcrit
        fit$DoF <- ffit$mdf
        if(CV){	
            fit$kappa <- ffit$kappa	
        }else{
             fit$kappa <- kappa
         }
    }
    fit$na.action <- "na.fail"
    # }}}
    if (method=="weib") fit$method <- "Weib" else fit$method <- "Splines"
    class(fit) <- "idm"
    fit$runtime <- proc.time()-ptm
    fit
}
