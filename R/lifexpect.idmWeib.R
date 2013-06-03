# 0 : health state
# 1 : illness state
# 2 : death state

##### Fonction qui calcule les esperances de vie avec leurs intervalles de confiance au temps s
lifexpect.idmWeib <- function(object,s,Z01,Z02,Z12,nsim=1000,CI=TRUE,...) {
	xx <- object
	if (inherits(xx,"idmWeib")) {
		nvar01 <- xx$NC[1]
		nvar02 <- xx$NC[2]
		nvar12 <- xx$NC[3]
		if(missing(Z01) && nvar01>0) warning("value(s) used for covariate(s) on transition 01: 0 \n")
		if(missing(Z02) && nvar02>0) warning("value(s) used for covariate(s) on transition 02: 0 \n")
		if(missing(Z12) && nvar12>0) warning("value(s) used for covariate(s) on transition 12: 0 \n")
		a01 <- xx$modelPar[1]
		b01 <- xx$modelPar[2]
		a02 <- xx$modelPar[3]
		b02 <- xx$modelPar[4]
		a12 <- xx$modelPar[5]
		b12 <- xx$modelPar[6]
		if (!missing(Z01))  {
			if (length(Z01) != nvar01) {
				stop("The length of the Z01 arguments must match the number of covariates on the ij transition.")
			}else{
				beta01 <- xx$coef[1:nvar01]
				bZ01 <- t(beta01)%*%Z01
			}
		}else{
			if (nvar01 != 0) { beta01 <- xx$coef[1:nvar01] }
			else { beta01 <- NULL }
			bZ01 <- 0
		}

		if (!missing(Z02))  {
			if (length(Z02) != nvar02) {
				stop("The length of the Z02 arguments must match the number of covariates on the ij transition.")
			}else{
				beta02 <- xx$coef[(nvar01+1):(nvar01+nvar02)]
				bZ02 <- t(beta02)%*%Z02
			}
		}else{
			if (nvar02 != 0) { beta02 <- xx$coef[(nvar01+1):(nvar01+nvar02)] }
			else { beta02 <- NULL }
			bZ02 <- 0
		}

		if (!missing(Z12))  {
			if (length(Z12) != nvar12) {
				stop("The length of the Z12 arguments must match the number of covariates on the ij transition.")
			}else{
				beta12 <- xx$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12)]
				bZ12 <- t(beta12)%*%Z12
			}
		}else{
			if (nvar12 != 0) { beta12 <- xx$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12)] }
			else { beta12 <- NULL }
			bZ12 <- 0
		}
	
		res <- lifexpect0.idmWeib(s,a01,1/b01,a02,1/b02,a12,1/b12,bZ01,bZ02,bZ12)
		
		if (CI==TRUE) {
		### CI prediction by Monte-Carlo
			Vmean <- c(sqrt(a01),sqrt(b01),sqrt(a02),sqrt(b02),sqrt(a12),sqrt(b12),beta01,beta02,beta12) # vector of parameters
			Mcov = xx$V
			# une simulation pour chaque element d'une liste
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
			X <- rmvnorm(nsim,Vmean,Mcov) 
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

		Xres1 <- mapply(function(Xa01,Xb01,Xa02,Xb02,Xa12,Xb12,XbZ01,XbZ02,XbZ12) lifexpect0.idmWeib(s,a01=Xa01,b01=Xb01,a02=Xa02,b02=Xb02,a12=Xa12,b12=Xb12,bZ01=XbZ01,bZ02=XbZ02,bZ12=XbZ12),Xa01,Xb01,Xa02,Xb02,Xa12,Xb12,XbZ01,XbZ02,XbZ12)
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
		return(list(HLE=c(res$HLE,Xres4[1,]),
		LE0=c(res$LE0,Xres4[2,]),
		LE1=c(res$LE1,Xres4[3,])))	
		}
		else 
			return(res)
  		
	}

}

lifexpect0.idmWeib <- function(s,a01,b01,a02,b02,a12,b12,bZ01=0,bZ02=0,bZ12=0) {
	ET12 = integrate(
		f=function(x) {
               	S.weib(s,x,a12,b12,bZ12)
               	},s,Inf)
	ET0. = integrate(f=function(x) {
		S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)
		},s,Inf)
	ET01 = integrate(f=function(x) {
		sapply(x,function(x) {integrate(f=function(y)
		{
		 S.weib(s,y,a01,b01,bZ01)*S.weib(s,y,a02,b02,bZ02)*iweibull(y,a01,b01,bZ01)*S.weib(y,x,a12,b12,bZ12)
		}
		,lower=s,upper=x)$value})
               },s,Inf)
return(list(HLE=ET0.$value,LE0=ET01$value+ET0.$value,LE1=ET12$value))

}



