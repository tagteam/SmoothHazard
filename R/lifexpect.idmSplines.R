# 0 : health state
# 1 : illness state
# 2 : death state

##### Fonction qui calcule les esperances de vie avec leurs intervalles de confiance au temps s
lifexpect.idmSplines <- function(xx,s,Z01,Z02,Z12,nsim,CI,...) {
		nvar01 <- xx$NC[1]
		nvar02 <- xx$NC[2]
		nvar12 <- xx$NC[3]
		if(missing(Z01) && nvar01>0) warning("value(s) used for covariate(s) on transition 01: 0 \n")
		if(missing(Z02) && nvar02>0) warning("value(s) used for covariate(s) on transition 02: 0 \n")
		if(missing(Z12) && nvar12>0) warning("value(s) used for covariate(s) on transition 12: 0 \n")
		nz01 <- xx$nknots01
		nz02 <- xx$nknots02
		nz12 <- xx$nknots12
		zi01 <- xx$knots01
		zi02 <- xx$knots02
		zi12 <- xx$knots12
		the01 <- xx$theta01
		the02 <- xx$theta02
		the12 <- xx$theta12

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

		res <- lifexpect0(s,zi01,nz01,the01^2,zi12,nz12,the12^2,zi02,nz02,the02^2,bZ01,bZ12,bZ02)

		if (CI==TRUE) {
		### CI prediction by Monte-Carlo
			Vmean <- c(the01,the02,the12,beta01,beta02,beta12) # vector of parameters
			Mcov = xx$V
			# une simulation pour chaque element d'une liste
			Xtheta01 <- as.list(NULL)
			Xtheta02 <- as.list(NULL)
			Xtheta12 <- as.list(NULL)
			XbZ01 <- as.list(NULL)
			XbZ02 <- as.list(NULL)
			XbZ12 <- as.list(NULL)

			set.seed(21)
			X <- rmvnorm(nsim,Vmean,Mcov) 
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

		Xres1 <- mapply(function(Xtheta01,Xtheta12,Xtheta02,XbZ01,XbZ12,XbZ02) lifexpect0(s,zi01,nz01,Xtheta01,zi12,nz12,Xtheta12,zi02,nz02,Xtheta02,XbZ01,XbZ12,XbZ02),Xtheta01,Xtheta12,Xtheta02,XbZ01,XbZ12,XbZ02)
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
		return(list(life.in.0.expectancy=c(res$life.in.0.expectancy,Xres4[1,]),
		life.expectancy.nondis=c(res$life.expectancy.nondem,Xres4[2,]),
		life.expectancy.dis=c(res$life.expectancy.dem,Xres4[3,])))	
		}
		else 
                  return(res)
  		

              }

lifexpect0 <- function(s,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01=0,bZ12=0,bZ02=0) {
	ET12 = integrate(f=function(x) {
               Predict0.idmPl(s,x,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01,bZ12,bZ02)$p11
               },s,zi12[nz12+6])
  	ET0. = integrate(f=function(x) {
               Predict0.idmPl(s,x,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01,bZ12,bZ02)$p00
               },s,zi02[nz02+6])
	ET01 = integrate(f=function(x) {
               Predict0.idmPl(s,x,zi01,nz01,the01,zi12,nz12,the12,zi02,nz02,the02,bZ01,bZ12,bZ02)$p01
               },s,zi01[nz01+6])
return(list(life.in.0.expectancy=ET0.$value,life.expectancy.nondis=ET01$value+ET0.$value,life.expectancy.dis=ET12$value))

}

