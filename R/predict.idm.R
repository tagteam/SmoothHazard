# 0 : health state
# 1 : illness state
# 2 : death state
#' Predictions for an illness-death model using either a penalized likelihood
#' approach or a Weibull parametrization.
#' 
#' Predict transition probabilities and cumulative probabilities from an object
#' of class \code{idmSplines} with confidence intervals are calculated.
#' 
#' @param object an \code{idm} class objects returned by a call to the
#' \code{\link{idm}} function
#' @param s time point at which prediction is made.
#' @param t time horizon for prediction.
#' @param newdata A data frame with covariate values for prediction.  
#' @param nsim number of simulations for the confidence intervals
#' calculations.  The default is 200.
#' @param seed Seed passed to \code{set.seed} for Monte Carlo
#' simulation of confidence intervals.
#' @param conf.int Logical: with (\code{TRUE}) or without
#' (\code{FALSE}) confidence intervals for the predicted values. The
#' default is \code{TRUE}.
#' @param level Value between 0 and 1, the level of confidence. Default is .95.
#' @param lifeExpect Logical. If \code{TRUE} compute life
#' expectancies, i.e., \code{t=Inf}.
#' @param ... other parameters.
#' @return a list containing the following predictions with pointwise
#' confidence intervals: \item{p00}{the transition probability \eqn{p_{00}}.}
#' \item{p01}{the transition probability \eqn{p_{01}}.} \item{p11}{the
#' transition probability \eqn{p_{11}}.} \item{p12}{the transition probability
#' \eqn{p_{12}}.} \item{p02_0}{the probability of direct transition from state
#' 0 to state 2.} \item{p02_1}{the probability of transition from state 0 to
#' state 2 via state 1.} \item{p02}{transition probability \eqn{p_{02}}. Note
#' that \code{p02}=\code{p_02_0}+\code{p02_1}.} \item{F01}{the lifetime risk of
#' disease. \code{F01}=\code{p01}+\code{p02_1}.} \item{F0.}{the probability of
#' exit from state 0. \code{F0.}=\code{p02_0}+\code{p01}+\code{p02_1}.}
#' @author R: Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> and Thomas Alexander Gerds <tag@@biostat.ku.dk> Fortran:
#' Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' @seealso \code{\link{idm}} 
#' @keywords methods
#' @examples
#'
#' \dontrun{
#' d=simulateIDM(n = 100)
#' fit <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2+X3,
#'                formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3,
#'                data=d,conf.int=FALSE)
#' predict(fit,s=0,t=80,nsim=3,conf.int=FALSE,lifeExpect=F)
#' 
#' 
#' data(Paq1000)
#' library(prodlim)
#' fit <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#' formula01=Hist(time=list(l,r),event=dementia)~certif,data=Paq1000)
#' 
#' predict(fit,s=70,t=80,newdata=data.frame(certif=1))
#' predict(fit,s=70,lifeExpect=TRUE,t=80,newdata=data.frame(certif=1))
#' 
#' fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#' 		formula01=Hist(time=list(l,r),event=dementia)~certif,
#'                 formula12=~1,
#'                 method="Splines",
#' 		data=Paq1000)
#' 
#' pred <- predict(fit.splines,s=70,t=80,Z01=c(1),Z02=c(1))
#' pred
#' 
#' }
#'
#' @export
predict.idm <- function(object,s,t,newdata,nsim=200,seed=21,conf.int=TRUE,level=.95,lifeExpect=FALSE,...) {
    if (lifeExpect==TRUE) t <- Inf
    if (any(s>t)) {stop("You must respect the condition 's<t' to calculate p(s,t)")}
    x <- object
    Mvar = x$V # covariance matrix
    if (conf.int == TRUE){
        stopifnot(0<level && level<1)
        if (nsim < 2) stop("Need at least two simulations to construct confidence limits.")
        XbZ01 <- vector(nsim,mode="list")
        XbZ02 <- vector(nsim,mode="list")
        XbZ12 <- vector(nsim,mode="list")
    }
    if (missing(t) && lifeExpect==FALSE) stop("Argument t is missing.")
    if (lifeExpect==TRUE) t <- Inf
    if (missing(s)) stop("Argument s is missing.")
    ## if (missing(t) || is.infinite(t)) lifeExpect <- TRUE
    # if covariates: cov=c(cov1,cov2,cov3,...)
    nvar01 <- x$NC[1]
    nvar02 <- x$NC[2]
    nvar12 <- x$NC[3]
    if (!missing(newdata)){
        if (NROW(newdata)>1) stop("Argument newdata has more than one row\n.Currently this function works only for one covariate constallation at a time.")
        Z01 <- as.matrix(model.frame(formula=update.formula(formula(x$terms$Formula01),NULL~.),data=newdata))
        Z02 <- as.matrix(model.frame(formula=update.formula(formula(x$terms$Formula02),NULL~.),data=newdata))
        Z12 <- as.matrix(model.frame(formula=update.formula(formula(x$terms$Formula12),NULL~.),data=newdata))
    }else{
         vars <- unique(c(x$Xnames01,x$Xnames02,x$Xnames12))
         newdata <- data.frame(matrix(0,ncol=length(vars)))
         names(newdata) <- vars
         Z01 <- rep(0,length(x$Xnames01))
         Z02 <- rep(0,length(x$Xnames02))
         Z12 <- rep(0,length(x$Xnames12))
     }
    beta01 <- x$coef[1:nvar01]
    bZ01 <- Z01 %*% beta01
    beta02 <- x$coef[(nvar01+1):(nvar01+nvar02)]
    bZ02 <- Z02 %*% beta02
    beta12 <- x$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12)]
    bZ12 <- Z12 %*% beta12
    ## Splines
    if (x$method=="Splines"){
        nz01 <- x$nknots01
        nz02 <- x$nknots02
        nz12 <- x$nknots12
        zi01 <- x$knots01
        zi02 <- x$knots02
        zi12 <- x$knots12
        the01 <- x$theta01
        the02 <- x$theta02
        the12 <- x$theta12
        if (conf.int == TRUE){
            ### conf.int prediction by Monte-Carlo
            Vmean <- c(the01,the02,the12,beta01,beta02,beta12) # vector of estimates
            set.seed(seed)
            X <- mvtnorm::rmvnorm(nsim,Vmean,Mvar) 
            # 1 set of simulated parameters for each element of the list
            Xtheta01=as.list(X[,1:(nz01+2)]^2)
            Xtheta02=as.list(X[,(nz01+3):(nz01+nz02+4)]^2)
            Xtheta12=as.list(X[,(nz01+nz02+5):(nz01+nz02+nz12+6)]^2)
            XbZ01=X[,(nz01+nz02+nz12+7):(nz01+nz02+nz12+6+nvar01)] %*% t(Z01)
            XbZ02=X[,(nz01+nz02+nz12+6+nvar01+1):(nz01+nz02+nz12+6+nvar01+nvar02)] %*% t(Z02)
            XbZ12=X[,(nz01+nz02+nz12+6+nvar01+nvar02+1):(nz01+nz02+nz12+6+nvar01+nvar02+nvar12)] %*% t(Z12)
            if (lifeExpect==TRUE){
                sim.param <- do.call("rbind",
                                     lapply(1:nsim,function(i){
                                                lifexpect0.idmPl(s,zi01[[i]],nz01[[i]],Xtheta01[[i]],zi12[[i]],nz12[[i]],Xtheta12[[i]],zi02[[i]],nz02[[i]],Xtheta02[[i]],XbZ01[[i]],XbZ12[[i]],XbZ02[[i]])
                                            }))
            }else{
                 sim.param <- do.call("rbind",
                                      lapply(1:nsim,function(i){
                                                 Predict0.idmPl(s,t,zi01[[i]],nz01[[i]],Xtheta01[[i]],zi12[[i]],nz12[[i]],Xtheta12[[i]],zi02[[i]],nz02[[i]],Xtheta02[[i]],XbZ01[[i]],XbZ12[[i]],XbZ02[[i]])
                                             }))
             }
            q.lower <- (1-level)/2
            q.upper <- 1-q.lower
            ci <- apply(sim.param,2,quantile,c(q.lower,q.upper))
        }
        if (lifeExpect==TRUE){
            transprob <- lifexpect0.idmPl(s,zi01,nz01,the01^2,zi12,nz12,the12^2,zi02,nz02,the02^2,bZ01,bZ12,bZ02)
        }else{
             transprob <- Predict0.idmPl(s,t,zi01,nz01,the01^2,zi12,nz12,the12^2,zi02,nz02,the02^2,bZ01,bZ12,bZ02)
         }
    } else {
          a01 <- x$modelPar[1]
          b01 <- x$modelPar[2]
          a02 <- x$modelPar[3]
          b02 <- x$modelPar[4]
          a12 <- x$modelPar[5]
          b12 <- x$modelPar[6]
          if (conf.int==TRUE) {
              ### conf.int prediction by Monte-Carlo
              # vector of parameter estimates
              Vmean <- c(sqrt(a01),sqrt(b01),sqrt(a02),sqrt(b02),sqrt(a12),sqrt(b12),beta01,beta02,beta12) 
              set.seed(seed)
              X <- mvtnorm::rmvnorm(nsim,Vmean,Mvar) 
              as.list(1/(X[,2]^2))
              # set of simulated parameters for each element of the list
              Xa01=as.list(X[,1]^2)
              Xb01=as.list(1/(X[,2]^2))
              Xa02=as.list(X[,3]^2)
              Xb02=as.list(1/(X[,4]^2))
              Xa12=as.list(X[,5]^2)
              Xb12=as.list(1/(X[,6]^2))
              XbZ01=X[,(7:(6+nvar01))] %*% t(Z01)
              XbZ02=X[,((6+nvar01+1):(6+nvar01+nvar02))] %*% t(Z02)
              XbZ12=X[,((6+nvar01+nvar02+1):(6+nvar01+nvar02+nvar12))] %*% t(Z12)
              if (lifeExpect==TRUE){
                  sim.param <- do.call("rbind",
                                       lapply(1:nsim,function(i){
                                                  lifexpect0.idmWeib(s,a01=Xa01[[i]],b01=Xb01[[i]],a02=Xa02[[i]],b02=Xb02[[i]],a12=Xa12[[i]],b12=Xb12[[i]],bZ01=XbZ01[[i]],bZ02=XbZ02[[i]],bZ12=XbZ12[[i]])
                                              }))
              }else{
                   sim.param <- do.call("rbind",
                                        lapply(1:nsim,function(i){
                                                   Predict0.idmWeib(s,t,a01=Xa01[[i]],b01=Xb01[[i]],a02=Xa02[[i]],b02=Xb02[[i]],a12=Xa12[[i]],b12=Xb12[[i]],bZ01=XbZ01[[i]],bZ02=XbZ02[[i]],bZ12=XbZ12[[i]])
                                               }))
               }
              q.lower <- (1-level)/2
              q.upper <- 1-q.lower
              ci <- apply(sim.param,2,quantile,c(q.lower,q.upper))
          }
          if (lifeExpect==TRUE){
              transprob <- lifexpect0.idmWeib(s,a01,1/b01,a02,1/b02,a12,1/b12,bZ01,bZ02,bZ12)
          }else{
               transprob <- Predict0.idmWeib(s,t,a01,1/b01,a02,1/b02,a12,1/b12,bZ01,bZ02,bZ12)
           }
      }
    if (conf.int==TRUE){
        transprob <- lapply(data.frame(rbind(transprob,ci)),
                            function(x){x
                                        names(x) <- c("Estimate",paste("Lower",round(100*level),sep="."),paste("Upper",round(100*level),sep="."))
                                        x})
    }
    out <- list(transprob=transprob)
    out <- c(out,list(newdata=newdata))
    out <- c(out,list(s=s,t=t,conf.int=conf.int))
    class(out) <- "predict.idm"
    out
}

print.predict.idm <- function(x,digits=2,...){
    lifeExpect <- is.infinite(x$t)
    cat("Predictions of an irreversible illness-death model with states (0,1,2).\n")
    cat("For covariate values:\n")
    print(x$newdata,row.names=FALSE)
    cat("\n")
    fmt <- paste0("%1.", digits[[1]], "f")
    nix <- lapply(names(x$transprob),function(ntp){
                           tp <- x$transprob[[ntp]]
                           if (length(tp)==3)
                               paste(sprintf(tp[[1]],fmt=fmt)," [",sprintf(tp[[2]],fmt=fmt),",",sprintf(tp[[3]],fmt=fmt),"]",sep="")
                           else
                               sprintf(tp[[1]],fmt=fmt)
                       })
    px <- data.frame(Parameter=names(x$transprob),Estimate=do.call(rbind,nix))
    rownames(px) <- NULL
    if (lifeExpect==TRUE){
        cat("Remaining life expected sojourn times (starting at time ",x$s,"):\n\n")
        print(cbind("State"=c(0,1,2),px[px$Parameter %in% c("LE.0","LE.nondiseased","LE.diseased"),]),row.names=FALSE)
    }else{
         cat("For a subject in state '0' at time ",x$s,",\npredicted state occupation probability at time ",x$t,":\n\n",sep="")
         print(cbind("State"=c(0,1,2),px[px$Parameter %in% c("p00","p01","p02"),]),row.names=FALSE)
         cat("\nFor a subject in state '0' at time ",x$s,",\npredicted probability of exit from state 0 until time ",x$t,":\n\n",sep="")
         print(cbind("Path"=c("via state 1","any"),px[px$Parameter %in% c("F01","F0."),]),row.names=FALSE)
         cat("\nFor a subject in state '1' at time ",x$s,",\npredicted state occupation probability at time ",x$t,":\n\n",sep="")
         print(cbind("State"=c(1,2),px[px$Parameter %in% c("p11","p12"),]),row.names=FALSE)
     }
    invisible(px)
}

Predict0.idmPl <- function(s,t,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01=0,bZ12=0,bZ02=0) {
    if (s>(min(zi01[nz01+6],zi02[nz02+6],zi12[nz12+6]))) {stop("argument s is off")}    
    if (any(t>zi12[nz12+6])) {stop("argument t is off")}
    p11 <- S.pl(s,t,zi12,nz12,the12,bZ12)
    p12 <- 1-p11
    p00 <- S.pl(s,t,zi01,nz01,the01,bZ01)*S.pl(s,t,zi02,nz02,the02,bZ02)
    p02_0 <- sapply(t,function(t) {integrate(f=function(x){S.pl(s,x,zi01,nz01,the01,bZ01)*S.pl(s,x,zi02,nz02,the02,bZ02)*susp(x,zi02,nz02,the02,bZ02)$intensity},lower=s,upper=t)$value})
    p01 <- sapply(t,function(t) {integrate(f=function(x){S.pl(s,x,zi01,nz01,the01,bZ01)*S.pl(s,x,zi02,nz02,the02,bZ02)*susp(x,zi01,nz01,the01,bZ01)$intensity*S.pl(x,t,zi12,nz12,the12,bZ12)},lower=s,upper=t)$value})
    p02_1 <- 1-p00-p02_0-p01
    p02 <- p02_0+p02_1
    c(p00=p00,p01=p01,p11=p11,p12=p12,p02_0=p02_0,p02_1=p02_1,p02=p02,F01=p01+p02_1,F0.=p02_0+p01+p02_1)
}

Predict0.idmWeib <- function(s,t,a01,b01,a02,b02,a12,b12,bZ01=0,bZ02=0,bZ12=0) {
    p11 = S.weib(s,t,a12,b12,bZ12)
    p12 = 1-p11
    p00 = S.weib(s,t,a01,b01,bZ01)*S.weib(s,t,a02,b02,bZ02)
    p02_0 = sapply(t,function(t) {integrate(f=function(x){S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a02,b02,bZ02)},lower=s,upper=t)$value })
    p01 = sapply(t,function(t) {integrate(f=function(x){S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a01,b01,bZ01)*S.weib(x,t,a12,b12,bZ12)},lower=s,upper=t)$value})
    p02_1 = 1-p00-p02_0-p01
    p02 = p02_0+p02_1
    ## return(list(p00=p00,p01=p01,p11=p11,p12=p12,p02_0=p02_0,p02_1=p02_1,p02=p02,F01=p01+p02_1,F0.=p02_0+p01+p02_1))
    c(p00=p00,
      p01=p01,
      p11=p11,
      p12=p12,
      p02_0=p02_0,
      p02_1=p02_1,
      p02=p02,
      F01=p01+p02_1,
      F0.=p02_0+p01+p02_1)
}

# a = shape parameter
# b = scale parameter
iweibull <- function(x,a,b,bZ=0) {
    res = (a/b) * (x/b)**(a-1) * exp(bZ)
    return(res) 
}

# S(s,t) = S(t)/S(s)
# if S(s)=0, S(s,t)=0
S.weib <- function(s,t,a,b,bZ=0) {	
    res <- 0
    St <- (1-pweibull(t,shape=a,scale=b))^(exp(bZ))
    Ss <- (1-pweibull(s,shape=a,scale=b))^(exp(bZ))
    if (length(s)==1){
        if (Ss==0){res <- 0}
        else{res <- St/Ss}
    }else{
         idx0 <- which(Ss==0)
         idx <- which(Ss!=0)
         res[idx0] <- 0 
         res[idx] <- St/Ss[idx]
     }
    return(res)
}

susp <- function(x,zi,nz,the,bZ=0) {
    gl=rep(0,length(x))   # risque cumule
    lam=rep(0,length(x))  # risque
    su=rep(0,length(x))   # survie
    TF=rep(0,length(x)) # T si z[i-1]<=x[.]<z[i], F sinon
    som=0
    for (i in 5:(nz+3)) {
        TF = ( (zi[i-1]<=x) & (x<zi[i]) )
        if (sum(TF) != 0) { 
            ind = which(TF) 
            mm3=rep(0,length(ind))
            mm2=rep(0,length(ind))
            mm1=rep(0,length(ind))
            mm=rep(0,length(ind))
            im3=rep(0,length(ind))
            im2=rep(0,length(ind))
            im1=rep(0,length(ind))
            im=rep(0,length(ind))
            j = i-1
            if (j>4) { 
                som = sum(the[1:(j-4)])
            }
            ht = x[ind]-zi[j] #
            htm = x[ind]-zi[j-1] #
            h2t = x[ind]-zi[j+2] #
            ht2 = zi[j+1]-x[ind] #
            ht3 = zi[j+3]-x[ind] #
            hht = x[ind]-zi[j-2] #
            h = zi[j+1]-zi[j]
            hh = zi[j+1]-zi[j-1]
            h2 = zi[j+2]-zi[j]
            h3 = zi[j+3]-zi[j]
            h4 = zi[j+4]-zi[j]
            h3m = zi[j+3]-zi[j-1]
            h2n = zi[j+2]-zi[j-1]
            hn= zi[j+1]-zi[j-2]
            hh3 = zi[j+1]-zi[j-3]
            hh2 = zi[j+2]-zi[j-2]
            mm3[ind] = ((4*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2[ind] = ((4*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4*h2t*htm*ht2)/(hh2*h2n*hh*h))+((4*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1[ind] = (4*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4*htm*ht*h2t)/(h3m*h2*h*h2n))+((4*ht3*ht*ht)/(h3m*h3*h2*h))
            mm[ind] = 4*(ht*ht*ht)/(h4*h3*h2*h)
            im3[ind] = (0.25*(x[ind]-zi[j-3])*mm3[ind])+(0.25*hh2*mm2[ind])+(0.25*h3m*mm1[ind])+(0.25*h4*mm[ind])
            im2[ind] = (0.25*hht*mm2[ind])+(h3m*mm1[ind]*0.25)+(h4*mm[ind]*0.25)
            im1[ind] = (htm*mm1[ind]*0.25)+(h4*mm[ind]*0.25)
            im[ind] = ht*mm[ind]*0.25
            gl[ind] = som +(the[j-3]*im3[ind])+(the[j-2]*im2[ind])+(the[j-1]*im1[ind])+(the[j]*im[ind])
            lam[ind] = (the[j-3]*mm3[ind])+(the[j-2]*mm2[ind])+(the[j-1]*mm1[ind])+(the[j]*mm[ind])
        } # fin if (sum(TF) != 0)
    } # fin for
    TF = (x>=zi[nz+3])
    if (sum(TF) != 0) {
        ind = which(TF)
        som = sum(the[1:(nz+2)])
        gl[ind] = som
        lam[ind] = 4*the[nz+2]/(zi[nz+3]-zi[nz+2])
    }
    TF = (x<zi[4])
    if (sum(TF) != 0) {
        ind = which(TF)
        gl[ind] = 0
        lam[ind] = 0
    }
    e = exp(bZ)
    lam=lam*e
    gl=gl*e
    su = exp(-gl)
    return(list(intensity=lam,cumul.intensity=gl,survival=su))
}

A <- function(s,t,zi,nz,the,bZ=0) {
    res=rep(0,length(t))
    TF = (t>=zi[length(zi)])
    ind = which(TF)
    if (sum(TF)!=0) {res[ind]=susp((zi[nz+6]-10^-5),zi,nz,the,bZ)$cumul.intensity-susp(s,zi,nz,the,bZ)$cumul.intensity}
    TF = (t<zi[length(zi)])
    ind = which(TF)
    if (sum(TF)!=0) {res[ind]=susp(t[ind],zi,nz,the,bZ)$cumul.intensity-susp(s,zi,nz,the,bZ)$cumul.intensity}
    return(res)
}

### Survival function with two time s, t
# S(s,t) = S(t)/S(s)
#        = exp(-A(s,t))
S.pl <- function(s,t,zi,nz,the,bZ=0) {
	if (length(t)>=length(s)){
		res=rep(0,length(t))
		TF = (t>zi[length(zi)])
		ind = which(TF)
		if (sum(TF)!=0) {res[ind]=0}
		TF = (t<=zi[length(zi)])
		ind = which(TF)
		if (sum(TF)!=0) {res[ind]=susp(t[ind],zi,nz,the,bZ)$survival/susp(s,zi,nz,the,bZ)$survival}
	}else{		
		res=rep(0,length(s))
		if (t>zi[length(zi)]) {res=0}
		else {res=susp(t,zi,nz,the,bZ)$survival/susp(s,zi,nz,the,bZ)$survival}
	}
	return(res)
}



lifexpect0.idmWeib <- function(s,a01,b01,a02,b02,a12,b12,bZ01=0,bZ02=0,bZ12=0) {
    ## print("lifexpect0.idmWeib")
    ET12 = integrate(
        f=function(x) {
            S.weib(s,x,a12,b12,bZ12)
        },s,Inf)
    ET0dot = integrate(f=function(x) {
                           S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)
                       },s,Inf)
    ET01 = integrate(f=function(x){
                         sapply(x,function(x){
                                    integrate(f=function(y){
                                                  S.weib(s,y,a01,b01,bZ01)*S.weib(s,y,a02,b02,bZ02)*iweibull(y,a01,b01,bZ01)*S.weib(y,x,a12,b12,bZ12)},
                                              lower=s,
                                              upper=x)$value})},s,Inf)
    c(LE.0=ET0dot$value,
      LE.nondiseased=ET01$value+ET0dot$value,
      LE.diseased=ET12$value)

}

lifexpect0.idmPl <- function(s,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01=0,bZ12=0,bZ02=0) {
    ET12 = integrate(f=function(x) {
                         Predict0.idmPl(s,x,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01,bZ12,bZ02)$p11
                     },s,zi12[nz12+6])
    ET0dot = integrate(f=function(x) {
                           Predict0.idmPl(s,x,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01,bZ12,bZ02)$p00
                       },s,zi02[nz02+6])
    ET01 = integrate(f=function(x) {
                         Predict0.idmPl(s,x,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01,bZ12,bZ02)$p01
                     },s,zi01[nz01+6])
    c(LE.0=ET0dot$value,
      LE.nondiseased=ET01$value+ET0dot$value,
      LE.diseased=ET12$value)
}
