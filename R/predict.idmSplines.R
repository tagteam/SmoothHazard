# 0 : health state
# 1 : illness state
# 2 : death state

##### Fonction qui calcule les predictions avec leurs intervalles de confiance entre les temps s et t
predict.idmSplines <- function(object,s,t,Z01,Z02,Z12,nsim=2000,CI=TRUE,...) {

# if covariates: cov=c(cov1,cov2,cov3,...)

	x <- object
	if (inherits(x,"idmSplines")) {
		nvar01 <- x$NC[1]
		nvar02 <- x$NC[2]
		nvar12 <- x$NC[3]
		if(missing(Z01) && nvar01>0) warning("value(s) used for covariate(s) on transition 01: 0 \n")
		if(missing(Z02) && nvar02>0) warning("value(s) used for covariate(s) on transition 02: 0 \n")
		if(missing(Z12) && nvar12>0) warning("value(s) used for covariate(s) on transition 12: 0 \n")
		nz01 <- x$nknots01
		nz02 <- x$nknots02
		nz12 <- x$nknots12
		zi01 <- x$knots01
		zi02 <- x$knots02
		zi12 <- x$knots12
		the01 <- x$theta01
		the02 <- x$theta02
		the12 <- x$theta12
		if (!missing(Z01))  {
			if (length(Z01) != x$NC[1]) {
				stop("The length of the Z01 arguments must match the number of covariates on the ij transition.")
			}else{
				beta01 <- x$coef[1:nvar01]
				bZ01 <- t(beta01)%*%Z01
			}
		}else{
			if (nvar01 != 0) { beta01 <- x$coef[1:nvar01] }
			else { beta01 <- NULL }
			bZ01 <- 0
		}

		if (!missing(Z02))  {
			if (length(Z02) != x$NC[2]) {
				stop("The length of the Z02 arguments must match the number of covariates on the ij transition.")
			}else{
				beta02 <- x$coef[(nvar01+1):(nvar01+nvar02)]
				bZ02 <- t(beta02)%*%Z02
			}
		}else{
			if (nvar02 != 0) { beta02 <- x$coef[(nvar01+1):(nvar01+nvar02)] }
			else { beta02 <- NULL }
			bZ02 <- 0
		}

		if (!missing(Z12))  {
			if (length(Z12) != x$NC[3]) {
				stop("The length of the Z12 arguments must match the number of covariates on the ij transition.")
			}else{
				beta12 <- x$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12)]
				bZ12 <- t(beta12)%*%Z12
			}
		}else{
			if (nvar12 != 0) { beta12 <- x$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12)] }
			else { beta12 <- NULL }
			bZ12 <- 0
		}
	
		res <- Predict0.idmPl(s,t,zi01,nz01,the01^2,zi12,nz12,the12^2,zi02,nz02,the02^2,bZ01,bZ12,bZ02)

		if (CI!=FALSE) {
		### CI prediction by Monte-Carlo
			Vmean <- c(the01,the02,the12,beta01,beta02,beta12) # vector of estimates
			Mvar = x$V # covariance matrix
			Xtheta01 <- as.list(NULL)
			Xtheta02 <- as.list(NULL)
			Xtheta12 <- as.list(NULL)
			XbZ01 <- as.list(NULL)
			XbZ02 <- as.list(NULL)
			XbZ12 <- as.list(NULL)
			
			set.seed(21)
			X <- rmvnorm(nsim,Vmean,Mvar) 
			# 1 set of simulated parameters for each element of the list
			
			for(i in (1:nsim) ) {
				Xtheta01[[i]]=X[i,1:(nz01+2)]^2
				Xtheta02[[i]]=X[i,(nz01+3):(nz01+nz02+4)]^2
				Xtheta12[[i]]=X[i,(nz01+nz02+5):(nz01+nz02+nz12+6)]^2
				
			}
			for(i in (1:nsim) ) {
				if (!missing(Z01)) {			
					XbZ01[[i]]=X[i,(nz01+nz02+nz12+7):(nz01+nz02+nz12+6+nvar01)] %*% Z01
				} else {
					XbZ01[[i]]=0
				}
				if (!missing(Z02)) {	
					XbZ02[[i]]=X[i,(nz01+nz02+nz12+6+nvar01+1):(nz01+nz02+nz12+6+nvar01+nvar02)] %*% Z02
				} else {
					XbZ02[[i]]=0
				}
				if (!missing(Z12)) {	
					XbZ12[[i]]=X[i,(nz01+nz02+nz12+6+nvar01+nvar02+1):(nz01+nz02+nz12+6+nvar01+nvar02+nvar12)] %*% Z12
				} else {
					XbZ12[[i]]=0
				}
			}
		Xres1 <- mapply(function(Xtheta01,Xtheta12,Xtheta02,XbZ01,XbZ12,XbZ02) Predict0.idmPl(s,t,zi01,nz01,Xtheta01,zi12,nz12,Xtheta12,zi02,nz02,Xtheta02,XbZ01,XbZ12,XbZ02),Xtheta01,Xtheta12,Xtheta02,XbZ01,XbZ12,XbZ02)
		Xres2 <- apply(Xres1,2,as.numeric) # transforme le type des elts de la matrice de 'list' en 'numeric'
		Xres3 <- t(apply(Xres2,1,sort)) # classement des valeurs
		iinf=((nsim/100)*2.5) + 1
		isup=(nsim/100)*97.5
		if (iinf != round(iinf)){
			delta <- (ceiling(iinf)-iinf < iinf-floor(iinf))
			iinf <- ceiling(iinf)*delta + floor(iinf)*(1-delta)
		}
		if (isup != round(isup)){
			delta <- (ceiling(isup)-isup < isup-floor(isup))
			isup <- ceiling(isup)*delta + floor(isup)*(1-delta)
		}
		Xres4 <- cbind(Xres3[,iinf],Xres3[,isup]) # 1ere colonne = bornes inf pour chaque valeur ; 2eme colonne = borne sup pour chaque valeur
		return(list(p00=c(res$p00,Xres4[1,]),p01=c(res$p01,Xres4[2,]),p11=c(res$p11,Xres4[3,]),
		p12=c(res$p12,Xres4[4,]),p02_0=c(res$p02_0,Xres4[5,]),p02_1=c(res$p02_1,Xres4[6,]),
		p02=c(res$p02,Xres4[7,]),F01=c(res$F01,Xres4[8,]),F0.=c(res$F0.,Xres4[9,])))	
		}
		else 
			return(res)	
	}
}


Predict0.idmPl <- function(s,t,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01=0,bZ12=0,bZ02=0) {
	p11 <- rep(0,length(t))
	p12 <- rep(0,length(t))
	p00 <- rep(0,length(t))
	p02_0 <- rep(0,length(t))
	p01 <- rep(0,length(t))
	p02_1 <- rep(0,length(t))
	p02 <- rep(0,length(t))
	if (all(s>t)) {stop("You must respect the condition 's<t' to calculate p(s,t)")} 
	if (s>(min(zi01[nz01+6],zi02[nz02+6],zi12[nz12+6]))) {stop("argument s is out of bornes")}    
	if (any(t>zi12[nz12+6])) {stop("argument t is out of bornes")}
	  	
			p11 <- S.pl(s,t,zi12,nz12,the12,bZ12)
			p12 <- 1-p11
			p00 <- S.pl(s,t,zi01,nz01,the01,bZ01)*S.pl(s,t,zi02,nz02,the02,bZ02)
			p02_0 <- sapply(t,function(t) {integrate(f=function(x)
			  {S.pl(s,x,zi01,nz01,the01,bZ01)*S.pl(s,x,zi02,nz02,the02,bZ02)*susp(x,zi02,nz02,the02,bZ02)$intensity}
				,lower=s,upper=t)$value})
			p01 <- sapply(t,function(t) {integrate(f=function(x)
			  {S.pl(s,x,zi01,nz01,the01,bZ01)*S.pl(s,x,zi02,nz02,the02,bZ02)*susp(x,zi01,nz01,the01,bZ01)$intensity*S.pl(x,t,zi12,nz12,the12,bZ12)}
			,lower=s,upper=t)$value})
			p02_1 <- 1-p00-p02_0-p01
			p02 <- p02_0+p02_1

return(list(p00=p00,p01=p01,p11=p11,p12=p12,p02_0=p02_0,p02_1=p02_1,p02=p02,F01=p01+p02_1,F0.=p02_0+p01+p02_1))
}



