### intensity.function.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb  6 2016 (08:47) 
## Version: 
## last-updated: Feb  6 2016 (09:33) 
##           By: Thomas Alexander Gerds
##     Update #: 14
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##'  M-spline estimate of the transition intensity function
##' and the cumulative transition intensity function
##' for survival and illness-death models
##'
##' The estimate of the transition intensity function is a linear
##' combination of M-splines and the estimate of the cumulative transition
##' intensity function is a linear combination of I-splines (the integral of a
##' M-spline is called I-spline). The coefficients \code{theta} are the same for
##' the M-splines and I-splines.
##' 
##' Important: the theta parameters returned by \code{idm} are in fact the square root of the splines coefficients.
##' See examples.
##' 
##' @title M-spline estimate of the transition intensity function
##' @param times Time points at which to estimate the intensity function
##' @param knots Knots for the M-spline
##' @param number.knots Number of knots for the M-splines (and I-splines see details)
##' @param theta The coefficients for the linear combination of M-splines (and I-splines see details)
##' @param linear.predictor
##' @return
##' - intensity : the transition intensity function evaluated at \code{times}
##' - cumulative.intensity: the cumulative transition intensity function evaluated at \code{times}
##' - survival: the survival function, i.e., exp(-cumulative.intensity)
##' 
##' @seealso \code{\link{shr}}, \code{\link{idm}} 
##' @examples 
##' @export 
#' @author R: Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> and Thomas Alexander Gerds <tag@@biostat.ku.dk> 
intensity.function <- function(times,knots,number.knots,theta,linear.predictor=0) {
    cumulative.intensity=rep(0,length(times))   # risque cumule
    intensity=rep(0,length(times))  # risque
    survival=rep(0,length(times))   # survie
    TF=rep(0,length(times)) # T si z[i-1]<=times[.]<z[i], F sinon
    som=0
    for (i in 5:(number.knots+3)) {
        TF = ( (knots[i-1]<=times) & (times<knots[i]) )
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
                som = sum(theta[1:(j-4)])
            }
            ht = times[ind]-knots[j] #
            htm = times[ind]-knots[j-1] #
            h2t = times[ind]-knots[j+2] #
            ht2 = knots[j+1]-times[ind] #
            ht3 = knots[j+3]-times[ind] #
            hht = times[ind]-knots[j-2] #
            h = knots[j+1]-knots[j]
            hh = knots[j+1]-knots[j-1]
            h2 = knots[j+2]-knots[j]
            h3 = knots[j+3]-knots[j]
            h4 = knots[j+4]-knots[j]
            h3m = knots[j+3]-knots[j-1]
            h2n = knots[j+2]-knots[j-1]
            hn= knots[j+1]-knots[j-2]
            hh3 = knots[j+1]-knots[j-3]
            hh2 = knots[j+2]-knots[j-2]
            mm3[ind] = ((4*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2[ind] = ((4*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4*h2t*htm*ht2)/(hh2*h2n*hh*h))+((4*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1[ind] = (4*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4*htm*ht*h2t)/(h3m*h2*h*h2n))+((4*ht3*ht*ht)/(h3m*h3*h2*h))
            mm[ind] = 4*(ht*ht*ht)/(h4*h3*h2*h)
            im3[ind] = (0.25*(times[ind]-knots[j-3])*mm3[ind])+(0.25*hh2*mm2[ind])+(0.25*h3m*mm1[ind])+(0.25*h4*mm[ind])
            im2[ind] = (0.25*hht*mm2[ind])+(h3m*mm1[ind]*0.25)+(h4*mm[ind]*0.25)
            im1[ind] = (htm*mm1[ind]*0.25)+(h4*mm[ind]*0.25)
            im[ind] = ht*mm[ind]*0.25
            cumulative.intensity[ind] = som +(theta[j-3]*im3[ind])+(theta[j-2]*im2[ind])+(theta[j-1]*im1[ind])+(theta[j]*im[ind])
            intensity[ind] = (theta[j-3]*mm3[ind])+(theta[j-2]*mm2[ind])+(theta[j-1]*mm1[ind])+(theta[j]*mm[ind])
        } # fin if (sum(TF) != 0)
    } # fin for
    TF = (times>=knots[number.knots+3])
    if (sum(TF) != 0) {
        ind = which(TF)
        som = sum(theta[1:(number.knots+2)])
        cumulative.intensity[ind] = som
        intensity[ind] = 4*theta[number.knots+2]/(knots[number.knots+3]-knots[number.knots+2])
    }
    TF = (times<knots[4])
    if (sum(TF) != 0) {
        ind = which(TF)
        cumulative.intensity[ind] = 0
        intensity[ind] = 0
    }
    e = exp(linear.predictor)
    intensity=intensity*e
    cumulative.intensity=cumulative.intensity*e
    survival = exp(-cumulative.intensity)
    return(list(intensity=intensity,cumulative.intensity=cumulative.intensity,survival=survival))
}

#----------------------------------------------------------------------
### intensity.function.R ends here
