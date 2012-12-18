summary.survPl <- function(object,conf.int=.95,digits=4,pvalDigits=4,eps=.0001, ...){

	if (!inherits(object,"survPl")) stop("Object must be of class 'survPl'")
	x <- object
#	if (x$istop == 1){
	if (x$converged == 1){
		cat("Survival model using a Penalized Likelihood on the hazard function.\n")
		cat("\n")
		cat("number of subject: ", x$N,"\n")
		cat("number of covariates: ", x$NC,"\n")
		if(length(x$na.action))cat("observation deleted due to missing: ",length(x$na.action),"\n")

		if(x$NC>0){
			wald <- (x$coef/x$se)**2
			z <- abs(qnorm((1 + conf.int)/2))
	
			out <- data.frame("Hazard ratio"=format(round(exp(x$coef),digits)),
			"Standard error"=format(round(x$se,digits)),
			"CI.95"=paste("[",format(round(exp(x$coef - z * x$se),2)),";",format(round(exp(x$coef + z * x$se),2)),"]",sep=""),
			"P-value"=format.pval(1 - pchisq(wald, 1),digits=pvalDigits,eps=eps))
			rownames(out) <- names(x$coef)
			print(out,row.names=T)
		}
	}

}
