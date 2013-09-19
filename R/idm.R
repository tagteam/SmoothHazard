idm <- function(formula01,
                formula02,
                formula12,
                data,
                maxiter=200,
                eps=c(5,5,3),
                n.knots=c(7,7,7),
                knots="equidistant",
                CV=FALSE,
                kappa=c(800000,200000,50000),
                igraph=1,
                intensities="Weib",
                print.iter=FALSE,
                subset,
                na.action = na.fail){
    
  # {{{ check formula
  call <- match.call()
  ptm <- proc.time()
  if(missing(formula01))stop("The argument formula01 needs a formula.")
  if(missing(formula02))stop("The argument formula02 needs a formula.")	
  if(class(formula01)!="formula")stop("The argument formula01 must be a class 'formula'.")	
  if(class(formula02)!="formula")stop("The argument formula02 must be a class 'formula'.")		
  if(missing(formula02)) formula02 <- formula01
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
  # {{{ extract response
  responseTrans <- model.response(m01)
  
  responseAbs <- model.response(m02)
  # {{{ extract covariates
  ## formula01
  x01 <- model.matrix(formula01,data=m01)[, -1, drop = FALSE]
  NC01 <- NCOL(x01)
  Xnames01 <- paste(colnames(x01),"01",sep="_")

  ## formula02
  x02 <- model.matrix(formula02,data=m02)[, -1, drop = FALSE]
  NC02 <- NCOL(x02)
  ## Xnames02 <- colnames(x02)
  Xnames02 <- paste(colnames(x02),"02",sep="_")

  ## formula12
  x12 <- model.matrix(formula12,data=m12)[, -1, drop = FALSE]
  NC12 <- NCOL(x12)
  ## Xnames12 <- colnames(x12)
  Xnames12 <- paste(colnames(x12),"12",sep="_")
  # }}}
  # {{{ prepare censored event times 
  N <- length(responseAbs[,"time"])

  isIntervalCensored <- attr(responseTrans,"cens.type")=="intervalCensored"
  truncated <- nchar(attr(responseAbs,"entry.type"))>1
  
  if (truncated==0){
    entrytime <- as.double(NULL)
  }else{
    entrytime <- as.double(responseAbs[,"entry"])
  }
	
  if (isIntervalCensored){
    Ltime <- as.double(responseTrans[,"L",drop=TRUE])
    Rtime <- as.double(responseTrans[,"R",drop=TRUE])
  }else{
    Ltime <- as.double(responseTrans[,"time",drop=TRUE])
    Rtime <- as.double(responseTrans[,"time",drop=TRUE])
    # Rtime <- rep(0,N)
  }
  
  idm <- responseTrans[,"status"]==(as.integer(isIntervalCensored)+1)
  idd <- responseAbs[,"status"]==1
 
  if(!(intensities %in% c("Weib","Splines"))) stop("The intensities argument must be 'Weib' or 'Splines'")
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
  
  if (intensities == "Weib"){
    #	cat("------ Program Weibull ------ \n")
    size1 <- NC01 + NC02 + NC12
    size2 <- size1^2
    size_V <- size1 + 6
    ffit <- .Fortran("idmWeib",
                     ## input
                     as.double(entrytime),               #
                     as.double(Ltime),                   #l=
                     as.double(Rtime),                   #r=
                     as.double(responseAbs[,"time"]),    #d=
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
                     converged=as.integer(0),
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
                     as.integer(print.iter),
                     V_tot=as.double(matrix(0,nrow=size_V,ncol=size_V)),
                     package="SmoothHazard")
    
    ## ===Fortran delivers===
    ## Variable name| Explanation|Dimension|Storage mode|Remark
    ## loglik|log-likelihood without and with covariate|length 2|double|
    ## basepar|Weibull parameters|length 2|double|
    ## regpar|Regression coefficients|length P01+P02+P12|double| 0--->1 then 0--->2 then 1--->2
    ## v|covariance matrix|length (P01+P02+P12)*(P01+P02+P12)|double|
    ## converged|1=converged,2=iter > maxiter ,3=no|length 1|integer|
    ## cv | value of convergence criteria 1:likelihood, 2:parameter est, 3: gradient |length 3 | double |
    ## niter | number of iterations used to converge | lenght 1 | integer |
    ## t|time to plot alpha(t) |length 100|double|
    ## a01|transition intensity function for 0--->1|length 100|double|
    ## a01_l|lower confidence band for a01|length 100|double|
    ## a01_u|upper confidence band for a01|length 100|double|
    ## a02|transition intensity function for 0--->2|length 100|double|
    ## a02_l|lower confidence band for a02|length 100|double|
    ## a02_u|upper confidence band for a02|length 100|double|
    ## a12|transition intensity function for 1--->2|length 100|double|
    ## a12_l|lower confidence band for a12|length 100|double|
    ## a12_u|upper confidence band for a12|length 100|double|
}else{
    #  	cat("------ Program Splines ------ \n")
    if (length(entrytime)>0)
        alltimes <- sort(unique(Ltime, Rtime,entrytime,responseAbs[,"time"]))
    else
        alltimes <- sort(unique(Ltime, Rtime,responseAbs[,"time"]))
    amax <- max(alltimes)
    amin <- min(alltimes)
    if (knots=="equidistant"){
        nknots01 <- n.knots[1]
        nknots02 <- n.knots[2]
        nknots12 <- n.knots[3]
        knots01 <- seq(amin,amax,(amax-amin)/(nknots01-1))
        knots02 <- seq(amin,amax,(amax-amin)/(nknots02-1))
        knots12 <- seq(amin,amax,(amax-amin)/(nknots12-1))
    }
    else{
        if (knots=="quantiles"){
            nknots01 <- n.knots[1]
            nknots02 <- n.knots[2]
            nknots12 <- n.knots[3]
            knots01 <- quantile(alltimes,seq(0,1,1/(nknots01-1)))
            knots02 <- quantile(alltimes,seq(0,1,1/(nknots02-1)))
            knots12 <- quantile(alltimes,seq(0,1,1/(nknots12-1)))
        }
        else{## user specified knots
            knots01 <- knots[[1]]
            knots02 <- knots[[2]]
            knots12 <- knots[[3]]
            nknots01 <- length(knots01)
            nknots02 <- length(knots02)
            nknots12 <- length(knots12)
        }
    }
    ## print(knots01)
    ## print(knots02)
    ## print(knots12)
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
                     as.double(responseAbs[,"time"]),    #d=
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
                     converged=as.integer(0),
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
                     as.integer(igraph),
                     CVcrit=as.double(0),
                     mdf=as.double(0),
                     theta01=as.double(rep(0,(nknots01+2))),
                     theta02=as.double(rep(0,(nknots02+2))),
                     theta12=as.double(rep(0,(nknots12+2))),
                     as.integer(print.iter),
                     V_tot=as.double(matrix(0,nrow=size_V,ncol=size_V)),
                     package="SmoothHazard")
}

  if (ffit$converged == 4){
      warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")    
  }

  if (ffit$converged == 2){
    warning("Model did not converge. You could change the 'maxit' parameter")
  }

  if (ffit$converged == 3){
    warning("Fisher information matrix non-positive definite.")
  }

  fit <- NULL
  if(intensities=="Weib"){
    weibullParameter <- ffit$basepar
  }

  NC <- c(NC01,NC02,NC12)
   
  fit$call <- call
  fit$loglik <- ffit$loglik
  fit$cv <- ffit$cv
  fit$niter <- ffit$niter
  fit$converged <- ffit$converged
  if(intensities=="Weib"){
    fit$modelPar <- weibullParameter
  }
  fit$N <- N
  fit$events1 <- sum(idm)
  fit$events2 <- sum(idd)
  fit$NC <- NC
  fit$responseAbs <- responseAbs
  fit$responseTrans <- responseTrans
  if(intensities=="Splines"){
    fit$time <- matrix(ffit$t,ncol=3) 
  }else{
    fit$time <- ffit$t 
  }
  fit$intensity01 <- ffit$a01
  fit$lowerIntensity01 <- ffit$a01_l
  fit$upperIntensity01 <- ffit$a01_u
  
  fit$intensity02 <- ffit$a02
  fit$lowerIntensity02 <- ffit$a02_l
  fit$upperIntensity02 <- ffit$a02_u
  
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
  if(intensities=="Weib"){
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
  
  if(intensities=="Splines"){
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
  fit$na.action <- na.action	
  # }}}
  if(intensities=="Weib"){
    class(fit) <- "idmWeib"
  }else{
    class(fit) <- "idmSplines"
  }  
	
  cost<-proc.time()-ptm
  ## cat("The program took", round(cost[3],2), "seconds \n")
  fit
}
