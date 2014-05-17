#' Print method for \code{idmWeib} objects
#' 
#' Print a summary of a fitted illness-death model using the Weibull approach.
#' 
#' 
#' @param x Class \code{idm} object, i.e. the result of a call to the
#' \code{\link{idm}} function with \code{intensities}="Weib".
#' @param conf.int Confidence level.
#' @param digits Number of digits to print.
#' @param pvalDigits Number of digits to print for p-values.
#' @param eps Passed to \code{format.pval}.
#' @param \dots Not used.
#' @author Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> 
#' @seealso \code{\link{summary.idmWeib}}, \code{\link{plot.idmWeib}}
#' @keywords methods
#' @examples
#' 
#' \dontrun{
#' data(Paq1000)
#' fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=t0)~certif,
#' 		formula01=Hist(time=list(l,r),event=dementia)~certif,
#'                 formula12=~1,
#'                 method="Splines",
#' 		data=Paq1000)
#' print(fit.splines)
#' 
#' }
#' @S3method print idm
print.idm <- function(x,conf.int=.95,digits=4,pvalDigits=4,eps=0.0001,...){
    cl <- x$call
    cat("Call:\n")
    dput(cl)
    cat("\n")
    if ( (x$converged[1]== 1)&(x$converged[2] %in% c(0,1))  ){
        if(sum(x$NC)>0){
            wald <- (x$coef/x$se)**2
            z <- abs(qnorm((1 + conf.int)/2))
            tmp <- data.frame("coef"=format(round(x$coef,digits)),
                              "SE coef"=format(round(x$se,digits)),
                              "HR"=format(round(x$HR,digits)),
                              "CI"=paste("[",format(round(exp(x$coef - z * x$se),2)),";",format(round(exp(x$coef + z * x$se),2)),"]",sep=""),
                              "Wald"=format(wald,digits),
                              "p-value"=format.pval(1 - pchisq(wald, 1),digits=pvalDigits,eps=eps))
        }
        tmp1 <- matrix(x$loglik,nrow=1)
        dimnames(tmp1) <- list(paste(ifelse(object$method=="Splines","Penalized",""),
                                     "log-likelihood",
                                     c("Without covariates", "With covariates")))
        cat("Illness-death regression model using",
            ifelse(object$method=="Splines"," M-spline approximations","Weibull parametrization"),
            "\n of the baseline transition intensities.\n")
        cat("\n")
        cat("number of subjects: ", x$N,"\n")
        cat("number of events '0 ->1': ", x$events1,"\n")
        cat("number of events '0 ->2' or '0 -> 1 -> 2': ", x$events2,"\n")
        cat("number of covariates: ", x$NC,"\n")
        print(x$na.action)
        if (x$method=="Splines"){
            tmp2 <- data.frame("transition01"=c(x$nknots01,x$kappa[1]),"transition02"=c(x$nknots02,x$kappa[2]),"transition12"=c(x$nknots12,x$kappa[3]))
            rownames(tmp2) <- c("knots","kappa")
            cat("\n")
            if(x$CV){
                cat("Smoothing parameters estimated by cross validation:\n")
                print(tmp2,row.names=T)
                cat("Cross validation criterion:",x$CVcrit,"\n")
                cat("DoF: ", formatC(-x$DoF, format="f",digits=2),"\n")
            }else{
                cat("Smoothing parameters:\n")
                print(tmp2,row.names=T)
            }
        }else{
            cat("Parameters of the Weibull distributions: 'S(t) = exp(-(b*t)^a)'\n")
            tmp <- matrix(x$modelPar,nrow=2)
            dimnames(tmp) <- list(c("shape (a)","scale (b)"),
                                  c("transition 0 -> 1",
                                    "transition 0 -> 2",
                                    "transition 1 -> 2"))
            prmatrix(tmp)
        }
        cat("\n")
        if(sum(x$NC)>0){
            print(tmp,row.names=T)            
        }else{
            cat("No covariates in the model.\n")
        }
        cat("\n")
        prmatrix(tmp1)
        cat("\n")
        cat("----\nModel converged.\n")
        cat("number of iterations: ", x$niter,"\n")
        cat("convergence criteria: parameters=", signif(x$cv[1],2), "\n")
        cat("                    : likelihood=", signif(x$cv[2],2), "\n") 
        cat("                    : second derivatives=", signif(x$cv[3],2), "\n")
    }
    if( (sum(x$NC)>0)&(x$converged[1]==1)&(x$converged[2]!=1) ){
        warning("The model did converge without covariates but did not converge with covariates","\n")
        switch(as.character(x$converged[2]),
               "2"={ warning("Maximum number of iterations reached.",call.=FALSE)},
               "3"={ warning("Model did not converge.",call.=FALSE)})
        cat("Log-likelihood without covariates: ",x$loglik[1], "\n")
    }
    if( (sum(x$NC)>0)&(x$converged[1]!=1)&(x$converged[2]==1) ){
        cat("The model did converge with covariates but did not converge without covariates","\n")
        switch(as.character(x$converged[1]),
               "2"={ warning("Maximum number of iterations reached.",call.=FALSE)},
               "3"={ warning("Model did not converge.",call.=FALSE)})	
	
          wald <- (x$coef/x$se)**2
          z <- abs(qnorm((1 + conf.int)/2))
          Xnames <- NULL
          if(x$NC[1]>0) Xnames <- c(Xnames,paste(x$Xnames01,"_01",sep=""))
          if(x$NC[2]>0) Xnames <- c(Xnames,paste(x$Xnames02,"_02",sep=""))
          if(x$NC[3]>0) Xnames <- c(Xnames,paste(x$Xnames12,"_12",sep=""))
          tmp <- data.frame("coef"=format(round(x$coef,digits)),"SE coef"=format(round(x$se,digits)),
                            "HR"=format(round(x$HR,digits)),
                            "CI"=paste("[",format(round(exp(x$coef - z * x$se),2)),";",format(round(exp(x$coef + z * x$se),2)),"]",sep=""),
                            "Wald"=format(wald,digits),"p-value"=format.pval(1 - pchisq(wald, 1),digits=pvalDigits,eps=eps))
          rownames(tmp) <- Xnames
          tmp1 <- matrix(x$loglik,nrow=1)
          dimnames(tmp1) <- list("Log likelihood", c("Without cov", "With cov"))		
        
          cat("Illness-death model: Results of Weibull regression for the intensity functions.\n")
          cat("\n")
          cat("number of subjects: ", x$N,"\n")
          cat("number of events '0-->1': ", x$events1,"\n")
          cat("number of events '0-->2' or '0-->1-->2': ", x$events2,"\n")
          cat("number of covariates: ", x$NC,"\n")
          ## if(length(x$na.action))cat("number of deleted observations due to missing: ",length(x$na.action),"\n")
                
          cat("\n")
          print(tmp,row.names=T)
          cat("\n")

          cat("Warning: no convergence without covariates","\n")
          
          prmatrix(tmp1)
          cat("\n")
          cat("Parameters of the Weibull distribution: 'S(t) = exp(-(b*t)^a)'\n")
          tmp <- matrix(x$modelPar,nrow=2)
          dimnames(tmp) <- list(c("a","b"),c("alpha01", "alpha02", "alpha12"))
          prmatrix(tmp)
          cat("\n")
          cat("----\nModel converged.\n")
          cat("number of iterations: ", x$niter,"\n")
          cat("convergence criteria: parameters=", signif(x$cv[1],2), "\n")
          cat("                    : likelihood=", signif(x$cv[2],2), "\n") 
          cat("                    : second derivatives=", signif(x$cv[3],2), "\n")

    }
    if( ((x$converged[1]!=1)&(x$converged[2]!=1)) ){
        warning("The model did not converge.","\n")
        switch(as.character(x$converged[1]),
               "2"={ warning("Maximum number of iterations reached.",call.=FALSE)},
               "3"={ warning("Fisher information matrix non-positive definite.",call.=FALSE)})	
    }
}
