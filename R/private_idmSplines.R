

susp <- function(x,zi,nz,the,bZ=0) {


 gl=rep(0,length(x))   # risque cumule
 lam=rep(0,length(x))  # risque
 su=rep(0,length(x))   # survie
 TF=rep(0,length(x)) # T si z[i-1]<=x[.]<z[i], F sinon
 som=0
	for (i in 5:(nz+3)) {
    TF = ( (zi[i-1]<=x) & (x<zi[i]) )
    if (sum(TF) != 0) { 
      ind = which(TF) 
      mm3=rep(0,length(ind))
      mm2=rep(0,length(ind))
      mm1=rep(0,length(ind))
      mm=rep(0,length(ind))
      im3=rep(0,length(ind))
      im2=rep(0,length(ind))
      im1=rep(0,length(ind))
      im=rep(0,length(ind))
			j = i-1
			if (j>4) { 
        			som = sum(the[1:(j-4)])
			}
			ht = x[ind]-zi[j] #
			htm = x[ind]-zi[j-1] #
			h2t = x[ind]-zi[j+2] #
			ht2 = zi[j+1]-x[ind] #
			ht3 = zi[j+3]-x[ind] #
			hht = x[ind]-zi[j-2] #
			h = zi[j+1]-zi[j]
			hh = zi[j+1]-zi[j-1]
			h2 = zi[j+2]-zi[j]
			h3 = zi[j+3]-zi[j]
			h4 = zi[j+4]-zi[j]
			h3m = zi[j+3]-zi[j-1]
			h2n = zi[j+2]-zi[j-1]
			hn= zi[j+1]-zi[j-2]
			hh3 = zi[j+1]-zi[j-3]
			hh2 = zi[j+2]-zi[j-2]
			mm3[ind] = ((4*ht2*ht2*ht2)/(h*hh*hn*hh3))
			mm2[ind] = ((4*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4*h2t*htm*ht2)/(hh2*h2n*hh*h))+((4*h2t*h2t*ht)/(hh2*h2*h*h2n))
			mm1[ind] = (4*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4*htm*ht*h2t)/(h3m*h2*h*h2n))+((4*ht3*ht*ht)/(h3m*h3*h2*h))
			mm[ind] = 4*(ht*ht*ht)/(h4*h3*h2*h)
			im3[ind] = (0.25*(x[ind]-zi[j-3])*mm3[ind])+(0.25*hh2*mm2[ind])+(0.25*h3m*mm1[ind])+(0.25*h4*mm[ind])
			im2[ind] = (0.25*hht*mm2[ind])+(h3m*mm1[ind]*0.25)+(h4*mm[ind]*0.25)
			im1[ind] = (htm*mm1[ind]*0.25)+(h4*mm[ind]*0.25)
			im[ind] = ht*mm[ind]*0.25
			gl[ind] = som +(the[j-3]*im3[ind])+(the[j-2]*im2[ind])+(the[j-1]*im1[ind])+(the[j]*im[ind])
			lam[ind] = (the[j-3]*mm3[ind])+(the[j-2]*mm2[ind])+(the[j-1]*mm1[ind])+(the[j]*mm[ind])
    } # fin if (sum(TF) != 0)
	} # fin for
  TF = (x>=zi[nz+3])
  if (sum(TF) != 0) {
    ind = which(TF)
    som = sum(the[1:(nz+2)])
		gl[ind] = som
		lam[ind] = 4*the[nz+2]/(zi[nz+3]-zi[nz+2])
  }
  TF = (x<zi[4])
  if (sum(TF) != 0) {
    ind = which(TF)
		gl[ind] = 0
		lam[ind] = 0
  }
  
 	e = exp(bZ)
  	lam=lam*e
  	gl=gl*e
	su = exp(-gl)
	
	return(list(intensity=lam,cumul.intensity=gl,survival=su))
}

A <- function(s,t,zi,nz,the,bZ=0) {
	res=rep(0,length(t))
	TF = (t>=zi[length(zi)])
	ind = which(TF)
	if (sum(TF)!=0) {res[ind]=susp((zi[nz+6]-10^-5),zi,nz,the,bZ)$cumul.intensity-susp(s,zi,nz,the,bZ)$cumul.intensity}
	TF = (t<zi[length(zi)])
	ind = which(TF)
	if (sum(TF)!=0) {res[ind]=susp(t[ind],zi,nz,the,bZ)$cumul.intensity-susp(s,zi,nz,the,bZ)$cumul.intensity}
	return(res)
}

### fonction de survie entre 2 temps s et t
# S(s,t) = S(t)/S(s)
#        = exp(-A(s,t))
S.pl <- function(s,t,zi,nz,the,bZ=0) {
	if (length(t)>=length(s)){
		res=rep(0,length(t))
		TF = (t>zi[length(zi)])
		ind = which(TF)
		if (sum(TF)!=0) {res[ind]=0}
		TF = (t<=zi[length(zi)])
		ind = which(TF)
		if (sum(TF)!=0) {res[ind]=susp(t[ind],zi,nz,the,bZ)$survival/susp(s,zi,nz,the,bZ)$survival}
	}else{		
		res=rep(0,length(s))
		if (t>zi[length(zi)]) {res=0}
		else {res=susp(t,zi,nz,the,bZ)$survival/susp(s,zi,nz,the,bZ)$survival}
	}
	return(res)
}

