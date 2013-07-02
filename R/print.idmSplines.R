print.idmSplines <- function(x,conf.int=.95,digits=4,pvalDigits=4,eps=0.0001,...){
	cl <- x$call
	cat("Call:\n")
	dput(cl)
	cat("\n")

#	if (x$istop == 1){
	if (x$converged == 1){
		if(sum(x$NC)>0){
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
		}
		tmp1 <- matrix(x$loglik,nrow=1)
		dimnames(tmp1) <- list("Penalized log likelihood", c("Without cov", "With cov"))		
		
		cat("Illness-death model using a penalized likelihood approach with splines approximation for the intensity functions.\n")
		cat("\n")
		cat("number of subjects: ", x$N,"\n")
		cat("number of events '0-->1': ", x$events1,"\n")
		cat("number of events '0-->2' or '0-->1-->2': ", x$events2,"\n")
		cat("number of subjects: ", x$N,"\n")
		cat("number of covariates: ", x$NC,"\n")
		if(length(x$na.action))cat("number of deleted observations due to missing: ",length(x$na.action),"\n")

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
	}else{
		switch(as.character(x$converged),
		"2"={ warning("Maximum number of iterations reached.",call.=FALSE)},
		"3"={ warning("Model did not converge.",call.=FALSE)})	
	}
}
