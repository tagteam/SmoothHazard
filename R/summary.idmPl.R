summary.idmPl <- function(object,conf.int=.95,digits=4,pvalDigits=4,eps=.0001, ...){

	if (!inherits(object,"idmPl")) stop("Object must be of class 'idmPl'")
	x <- object
	if (x$converged == 1){
		cat("illness using a Penalized Likelihood on the hazard function.\n")
		cat("\n")
		cat("number of subjects: ", x$N,"\n")
		cat("number of events '0-->1': ", x$events1,"\n")
		cat("number of events '0-->2 or 0-->1-->2': ", x$events2,"\n")
		cat("number of covariates: ", x$NC,"\n")
		if(length(x$na.action))cat("observation deleted due to missing: ",length(x$na.action),"\n")

		if(sum(x$NC)>0){
			wald <- (x$coef/x$se)**2
			z <- abs(qnorm((1 + conf.int)/2))
	
			out <- data.frame("Hazard ratio"=format(round(exp(x$coef),digits)),
			"Standard error"=format(round(x$se,digits)),
			"CI.95"=paste("[",format(round(exp(x$coef - z * x$se),2)),";",format(round(exp(x$coef + z * x$se),2)),"]",sep=""),
			"P-value"=format.pval(1 - pchisq(wald, 1),digits=pvalDigits,eps=eps))
			Xnames <- NULL
			if(x$NC[1]>0) Xnames <- c(Xnames,paste(x$Xnames01,"_01",sep=""))
			if(x$NC[2]>0) Xnames <- c(Xnames,paste(x$Xnames02,"_02",sep=""))
			if(x$NC[3]>0) Xnames <- c(Xnames,paste(x$Xnames12,"_12",sep=""))
			rownames(out) <- Xnames
			print(out,row.names=T)
		}
	}

}
