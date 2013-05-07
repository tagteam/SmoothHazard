print.shrWeib <- function(x,conf.int=.95,digits=4,pvalDigits=4,eps=0.0001,...){
	cl <- x$call
	cat("Call:\n")
	dput(cl)
	cat("\n")

	if (x$converged == 1){
		if(x$NC >0){
			wald <- (x$coef/x$se)**2
			z <- abs(qnorm((1 + conf.int)/2))
	
			tmp <- data.frame("coef"=format(round(x$coef,digits)),"SE coef"=format(round(x$se,digits)),
			"HR"=format(round(x$HR,digits)),
			"CI"=paste("[",format(round(exp(x$coef - z * x$se),2)),";",format(round(exp(x$coef + z * x$se),2)),"]",sep=""),
			"Wald"=format(wald,digits),"P-value"=format.pval(1 - pchisq(wald, 1),digits=pvalDigits,eps=eps))
			rownames(tmp) <- names(x$coef)
		}

		tmp1 <- matrix(x$loglik,nrow=1)
		dimnames(tmp1) <- list("Log likelihood", c("Without cov", "With cov"))		
		
		cat("Survival model using a parametric approach with a Weibull distribution for the hazard function.\n")
		cat("\n")
		## cat("number of subjects: ", x$N,"\n")
		## cat("number of events: ", x$events,"\n")
                print(x$modelResponse)
		cat("number of covariates: ", x$NC,"\n")
		if(length(x$na.action))cat("number of deleted observations due to missing: ",length(x$na.action),"\n")

		if(x$NC >0){
			cat("\n")
			print(tmp,row.names=T)
		}
		cat("\n")
		prmatrix(tmp1)
		cat("\n")
		cat("Parameters of the Weibull distribution: 'S(t) = exp(-(b*t)^a)'\n")
		cat("                   a = ",x$modelPar[1]," "," b = ",x$modelPar[2],"\n")
		cat("\n")
		cat("----\nModel converged.\n")	
		cat("number of iterations: ", x$niter,"\n")
		cat("Convergence criteria: parameters=", signif(x$cv[1],2), "\n")
		cat("                    : likelihood=", signif(x$cv[2],2), "\n") 
		cat("                    : second derivatives=", signif(x$cv[3],2), "\n")
	}else{
		switch(as.character(x$converged),
		"2"={ warning("Maximum number of iterations reached.",call.=FALSE)},
		"3"={ warning("Model did not converge.",call.=FALSE)})
	}

}
