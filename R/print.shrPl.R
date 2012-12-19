print.shrPl <- function(x,conf.int=.95,digits=4,pvalDigits=4,eps=0.0001,...){
	if (!inherits(x,"shrPl")) stop("Object must be of class 'shrPl'")

	cl <- x$call
	cat("Call:\n")
	dput(cl)
	cat("\n")

#	if (x$istop == 1){
	if (x$converged == 1){
### CT	
#		if(x$noVar!=1){
		if(x$NC >0){
### CT	
			wald <- (x$coef/x$se)**2
			z <- abs(qnorm((1 + conf.int)/2))
	
			tmp <- data.frame("coef"=format(round(x$coef,digits)),"SE coef"=format(round(x$se,digits)),
			"HR"=format(round(x$HR,digits)),
			"CI"=paste("[",format(round(exp(x$coef - z * x$se),2)),";",format(round(exp(x$coef + z * x$se),2)),"]",sep=""),
			"Wald"=format(wald,digits),"P-value"=format.pval(1 - pchisq(wald, 1),digits=pvalDigits,eps=eps))
			rownames(tmp) <- names(x$coef)
		}

		tmp1 <- matrix(x$loglik,nr=1)
		dimnames(tmp1) <- list("Log likelihood", c("Without cov", "With cov"))		
		
		cat("Survival model using a penalized likelihood approach with splines approximation for the hazard function.\n")
		cat("\n")
		cat("number of subjects: ", x$N,"\n")
		cat("number of events: ", x$events,"\n")
		cat("number of covariates: ", x$NC,"\n")
		cat("number of nodes: ", x$nknots,"\n")
		if(length(x$na.action))cat("number of deleted observations due to missing: ",length(x$na.action),"\n")

		cat("\n")
		if(x$irec == 1){
			cat("Smoothing parameters estimated by Cross validation: ",x$kappa,"\n")
			cat("Cross validation criterion:",x$CVcrit,"\n")
			cat("DoF: ", formatC(-x$DoF, format="f",dig=2),"\n")
		}else{
			cat("Smoothing parameters: ",x$kappa,"\n")
		}
### CT	
#		if(x$noVar!=1){
		if(x$NC >0){
### CT	
			cat("\n")
			print(tmp,row.names=T)
		}
		cat("\n")
		prmatrix(tmp1)
		cat("\n")
		cat("----\nModel converged.\n")	
		cat("number of iterations: ", x$niter,"\n")
		cat("convergence criteria: parameters=", signif(x$cv[1],2), "\n")
		cat("                    : likelihood=", signif(x$cv[2],2), "\n") 
		cat("                    : second derivatives=", signif(x$cv[3],2), "\n")
	}else{
		switch(as.character(x$converged),
		"2"={ warning("Maximum number of iterations reached.",call.=FALSE)},
		"3"={ warning("Model did not converge.",call.=FALSE)})
	}

}
