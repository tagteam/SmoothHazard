# 0 : health state
# 1 : illness state
# 2 : death state

##### Fonction qui calcule les predictions avec leurs intervalles de confiance entre les temps s et t
predict.idmWeib <- function(object,s,t,Z01,Z02,Z12,nsim=2000,CI=TRUE,...) {

# if covariates: cov=c(cov1,cov2,cov3,...)

	x <- object
	if (inherits(x,"idmWeib")) {
		nvar01 <- x$NC[1]
		nvar02 <- x$NC[2]
		nvar12 <- x$NC[3]
		if(missing(Z01) && nvar01>0) warning("value(s) used for covariate(s) on transition 01: 0 \n")
		if(missing(Z02) && nvar02>0) warning("value(s) used for covariate(s) on transition 02: 0 \n")
		if(missing(Z12) && nvar12>0) warning("value(s) used for covariate(s) on transition 12: 0 \n")
		a01 <- x$modelPar[1]
		b01 <- x$modelPar[2]
		a02 <- x$modelPar[3]
		b02 <- x$modelPar[4]
		a12 <- x$modelPar[5]
		b12 <- x$modelPar[6]
		if (!missing(Z01))  {
			if (length(Z01) != x$NC[1]) {
				stop("The length of the Z01 arguments must match the number of covariates on the ij transition.")
			}else{
				beta01 <- x$coef[1:nvar01]
				bZ01 <- t(beta01)%*%Z01
			}
		}else{
			if (nvar01 != 0) { beta01 <- x$coef[1:nvar01] 
			}else{ beta01 <- NULL 
			}
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
			if (nvar02 != 0){ beta02 <- x$coef[(nvar01+1):(nvar01+nvar02)] 
			}else{ beta02 <- NULL
			}
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
			if (nvar12 != 0){ beta12 <- x$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12)] 
			}else{ beta12 <- NULL 
			}
			bZ12 <- 0
		}

		res <- Predict0.idmWeib(s,t,a01,1/b01,a02,1/b02,a12,1/b12,bZ01,bZ02,bZ12)

		if (CI!=FALSE) {
		### CI prediction by Monte-Carlo
			Vmean <- c(sqrt(a01),sqrt(b01),sqrt(a02),sqrt(b02),sqrt(a12),sqrt(b12),beta01,beta02,beta12) # vector of estimates
			Mvar = x$V # covariance matrix
			Xa01 <- as.list(NULL)
			Xb01 <- as.list(NULL)
			Xa02 <- as.list(NULL)
			Xb02 <- as.list(NULL)
			Xa12 <- as.list(NULL)
			Xb12 <- as.list(NULL)
			XbZ01 <- as.list(NULL)
			XbZ02 <- as.list(NULL)
			XbZ12 <- as.list(NULL)
			
			set.seed(21)
			X <- rmvnorm(nsim,Vmean,Mvar) 
			# 1 set of simulated parameters for each element of the list
			
			for(i in (1:nsim) ) {
				Xa01[[i]]=X[i,1]^2
				Xb01[[i]]=1/(X[i,2]^2)
				Xa02[[i]]=X[i,3]^2
				Xb02[[i]]=1/(X[i,4]^2)
				Xa12[[i]]=X[i,5]^2
				Xb12[[i]]=1/(X[i,6]^2)
			}
			for(i in (1:nsim) ) {
				if (!missing(Z01)) {			
					XbZ01[[i]]=X[i,(7:(6+nvar01))] %*% Z01
				} else {
					XbZ01[[i]]=0
				}
				if (!missing(Z02)) {	
					XbZ02[[i]]=X[i,((6+nvar01+1):(6+nvar01+nvar02))] %*% Z02
				} else {
					XbZ02[[i]]=0
				}
				if (!missing(Z12)) {	
					XbZ12[[i]]=X[i,((6+nvar01+nvar02+1):(6+nvar01+nvar02+nvar12))] %*% Z12
				} else {
					XbZ12[[i]]=0
				}
			}
		Xres1 <- mapply(function(Xa01,Xb01,Xa02,Xb02,Xa12,Xb12,XbZ01,XbZ02,XbZ12) Predict0.idmWeib(s,t,a01=Xa01,b01=Xb01,a02=Xa02,b02=Xb02,a12=Xa12,b12=Xb12,bZ01=XbZ01,bZ02=XbZ02,bZ12=XbZ12),Xa01,Xb01,Xa02,Xb02,Xa12,Xb12,XbZ01,XbZ02,XbZ12)
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

Predict0.idmWeib <- function(s,t,a01,b01,a02,b02,a12,b12,bZ01=0,bZ02=0,bZ12=0) {
	p11=rep(0,length(t))
	p12=rep(0,length(t))
	p00=rep(0,length(t))
	p02_0=rep(0,length(t))
	p01=rep(0,length(t))
	p02_1=rep(0,length(t))
	p02=rep(0,length(t))
 	if (all(s>t)) {stop("You must respect the condition 's<t' to calculate p(s,t)")
	}else{
		p11 = S.weib(s,t,a12,b12,bZ12)
		p12 = 1-p11
		p00 = S.weib(s,t,a01,b01,bZ01)*S.weib(s,t,a02,b02,bZ02)
		p02_0 = sapply(t,function(t) {integrate(f=function(x)
		{S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a02,b02,bZ02)}
		,lower=s,upper=t)$value })
		p01 = sapply(t,function(t) {integrate(f=function(x)
		{S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a01,b01,bZ01)*S.weib(x,t,a12,b12,bZ12)}
		,lower=s,upper=t)$value})
		p02_1 = 1-p00-p02_0-p01
		p02 = p02_0+p02_1
 	}
return(list(p00=p00,p01=p01,p11=p11,p12=p12,p02_0=p02_0,p02_1=p02_1,p02=p02,F01=p01+p02_1,F0.=p02_0+p01+p02_1))
}


### fonction d'intensité de transition
# a = shape parameter
# b = scale parameter
# bz = (vecteur des coeffs de régression)^T * (vecteur des variables) 
iweibull <- function(x,a,b,bZ=0) {
  res = (a/b) * (x/b)**(a-1) * exp(bZ)
  return(res) 
}

#############################################


### fonction de survie entre 2 temps s et t
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



