##' Function to generate an illness-death model for simulation.
##'
##' Based on the functionality of the lava PACKAGE the function generates
##' a latent variable model (latent illtime, waittime and lifetime)
##' and censoring mechanism (censtime, inspection1,inspection2,...,inspectionK).
##' 
##' The function \code{\link{sim.idmModel}} then simulates
##' right censored lifetimes and interval censored illness times.
##'
##' @title Generate illness-death model objects
##' @param scale.illtime Weilbull scale for latent illness time
##' @param shape.illtime Weilbull shape for latent illness time
##' @param scale.lifetime Weilbull scale for latent life time
##' @param shape.lifetime Weilbull shape for latent life time
##' @param scale.waittime Weilbull scale for latent life time
##' @param shape.waittime Weilbull shape for latent life time
##' @param scale.censtime Weilbull scale for censoring time
##' @param shape.censtime Weilbull shape for censoring time
##' @param K Number of intervals
##' @param schedule Mean of the waiting time between adjacent inspections.
##' @param punctuality Standard deviation of waiting time between inspections.
##' @examples
##' library(lava)
##' # generate illness-death model based on exponentially
##' # distributed times
##' m <- idmModel()
##' sim(m,6)
##'
##' # Estimate the parameters of the Weibull models
##' # based on the uncensored exact event times
##' # and the uncensored illstatus.
##' set.seed(18)
##' d <- sim(m,10,latent=TRUE)
##' d$uncensored.status <- 1
##' f <- idm(formula01=Hist(time=illtime,event=illstatus)~1,
##'     formula02=Hist(time=lifetime,event=uncensored.status)~1,data=d,conf.int=FALSE)
##' print(f)
##'
##' # Lower the rate of the 0->2 and 0->1 transitions
##' # increase the rate of the 1->2 transition
##' # and also lower the censoring rate
##' m <- idmModel(scale.lifetime=1/2000,scale.waittime=1/30,scale.illtime=1/1000,scale.censtime=1/1000)
##' set.seed(18)
##' d <- sim(m,100,latent=TRUE)
##' d$uncensored.status <- 1
##' 
##' f <- idm(formula01=Hist(time=observed.illtime,event=illstatus)~1,
##'     formula02=Hist(time=observed.lifetime,event=uncensored.status)~1,data=d,conf.int=FALSE)
##' print(f)
##'
##' # Estimate based on the right censored observations
##' fc <- idm(formula01=Hist(time=illtime,event=seen.ill)~1,
##'     formula02=Hist(time=observed.lifetime,event=seen.exit)~1,data=d,conf.int=FALSE)
##' print(fc)
##'
##' # Estimate based on interval censored and right censored observations
##' fi <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~1,
##'     formula02=Hist(time=observed.lifetime,event=seen.exit)~1,data=d,conf.int=FALSE)
##' print(fi)
##' 
##' # Estimation of covariate effects:
##' # diabetes, systolic bloodpressure, good cholesterol (hdl)
##' m <- idmModel(shape.waittime=2,scale.lifetime=1/2000,scale.waittime=1/30,scale.illtime=1/1000,scale.censtime=1/1000)
##' distribution(m,"diabetes") <- binomial.lvm(p=0.3)
##' distribution(m,"sbp") <- normal.lvm(mean=120,sd=20)
##' distribution(m,"hdl") <- normal.lvm(mean=50,sd=20)
##' regression(m,to="latent.illtime",from="diabetes") <- 1.7
##' regression(m,to="latent.illtime",from="sbp") <- 0.07
##' regression(m,to="latent.illtime",from="hdl") <- -0.1
##' regression(m,to="latent.waittime",from="diabetes") <- 1.8
##' regression(m,to="latent.lifetime",from="diabetes") <- 0.7
##' set.seed(21)
##' d <- sim(m,500,latent=TRUE)
##' head(d)
##'
##' # Estimation based on uncensored data
##' d$uncensored.status <- 1
##' # uncensored data
##' F1 <- idm(formula01=Hist(time=illtime,event=illstatus)~diabetes+sbp+hdl,
##'     formula02=Hist(time=lifetime,event=uncensored.status)~diabetes+sbp+hdl,
##'     data=d,conf.int=FALSE)
##' print(F1)
##'
##' # Estimation based on right censored data
##' F2 <- idm(formula01=Hist(time=illtime,event=seen.ill)~diabetes+sbp+hdl,
##'     formula02=Hist(time=observed.lifetime,event=seen.exit)~diabetes+sbp+hdl,
##'     data=d,conf.int=FALSE)
##' print(F2)
##' 
##' # Estimation based on interval censored and right censored data
##' F3 <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~diabetes+sbp+hdl,
##'     formula02=Hist(time=observed.lifetime,event=seen.exit)~diabetes+sbp+hdl,
##'     data=d,conf.int=FALSE)
##' print(F3)
##' cbind(uncensored=F1$coef,right.censore=F2$coef,interval.censored=F3$coef)
##' 
##' @return A latent variable model object \code{lvm}
##' @author Thomas Alexander Gerds
##' @export
idmModel <- function(scale.illtime=1/100,
                     shape.illtime=1,
                     scale.lifetime=1/100,
                     shape.lifetime=1,
                     scale.waittime=1/100,
                     shape.waittime=1,
                     scale.censtime=1/100,
                     shape.censtime=1,
                     K=5,
                     schedule=10,
                     punctuality=5){
    
    ## illness-death-model
    ##
    ## model waiting time in state 0
    ## based on latent "ill" and "life" times
    ## separate model for the waiting time in state "ill"
    idm <- lava::lvm()
    lava::latent(idm) <- c("latent.illtime","latent.lifetime","latent.waittime","censtime")
    lava::distribution(idm,"latent.illtime") <- lava::coxWeibull.lvm(scale=scale.illtime,
                                                                     shape=shape.illtime)
    lava::distribution(idm,"latent.lifetime") <- lava::coxWeibull.lvm(scale=scale.lifetime,
                                                                      shape=shape.lifetime)
    ## idm <- lava::eventTime(idm,illtime~min(latent.illtime=1,latent.lifetime=0),"ill")
    ## idm <- lava::eventTime(idm,lifetime~min(illtime=1,latent.lifetime=2),"event")
    lava::distribution(idm,"latent.waittime") <- lava::coxWeibull.lvm(scale=scale.waittime,
                                                                      shape=shape.waittime)
    ## transform(m,lifetime~latent.lifetime+latent.illtime+ill+latent.waittime) <- function(x){
    ## if (x$ill==1)
    ## x$latent.illtime+x$latent.waittime
    ## else
    ## x$latent.lifetime
    ## }
    ## ===============================================
    ## right censoring time
    ## ===============================================
    lava::distribution(idm,"censtime") <- lava::coxWeibull.lvm(shape=shape.censtime,
                                                               scale=scale.censtime)
    if (K>0)
        for (k in 1:K){
            lava::distribution(idm,paste("inspection",k,sep="")) <- lava::normal.lvm(mean=schedule,sd=punctuality)
        }
    ## }
    class(idm) <- c("idmModel","lvm")
    ## class(idm) <- c("lvm")
    idm
}
##' Function to simulate illness-death model data
##'
##' Based on the functionality of the lava PACKAGE 
##' @title Simulate illness-death model data
##' @param x An \code{idmModel} object as obtained with
##' \code{idmModel}
##' @param n Number of observations
##' @param illness.known.at.death
##' @param compliance Probability of missing an inspection time.
##' @param latent if TRUE keep the latent event times
##' @param keep.inspectiontimes if \code{TRUE} keep the inspection
##' times.
##' @param ... Extra arguments given to \code{sim}
##' @return A data set with interval censored observations from an illness-death model
##' @examples
##' example(idmModel)
##' help(idmModel)
##' @author Thomas Alexander Gerds
##' @importFrom lava sim
##' @S3method sim idmModel
sim.idmModel <- function(x,
                         n,
                         illness.known.at.death=TRUE,
                         compliance=1,
                         latent=FALSE,
                         keep.inspectiontimes=FALSE,
                         ...){
    # simulate latent data
    class(x) <- "lvm"
    dat <- lava::sim(x,n=n,...)
    # construct illtime and true illness status
    dat$illtime <- dat$latent.illtime
    dat$illstatus <- 1*(dat$illtime<=dat$latent.lifetime)
    dat$illtime[dat$illtime>dat$latent.lifetime] <- 0
    # construct lifetime
    # for ill subjects as the sum of the time to illness (illtime) and
    # the time spent in the illness state (waittime)
    dat$lifetime <- dat$latent.lifetime
    dat$lifetime[dat$illstatus==1] <- dat$illtime[dat$illstatus==1]+dat$latent.waittime[dat$illstatus==1]
    # interval censored illtime
    ipos <- grep("inspection[0-9]+",names(dat))
    if (length(ipos)>0) {
        # compute inspection times
        # make sure all inspection times are in the future
        # of the previous inspection time
        iframe <- dat[,ipos]
        iframe <- do.call("rbind",lapply(1:n,function(i){
            cumsum(pmax(unlist(iframe[i,]),0))
        }))
        interval <- do.call("rbind",lapply(1:n,function(i){
            ## remove duplicates
            itimes <- unique(iframe[i,])
            ## remove inspections where compliance is
            ## sampled as zero
            if (compliance<1 & compliance>0){
                comp <- rbinom(length(itimes),1,compliance)
                itimes <- itimes[comp==1]
            }
            ## remove inspection times that are 
            ## larger than the individual lifetime
            ltime <- dat$lifetime[i]
            itimes <- itimes[itimes<ltime]
            ## set censoring time to the maximum of the
            ## random censoring time and the largest
            ## inspection time which is smaller than
            ## the lifetime
            ctime <- max(dat$censtime[i],itimes[length(itimes)])
            ## and add the individual lifetime as the largest
            ## inspection time
            itimes <- c(itimes,ltime)
            ## find interval where illness happens.
            ## if illness happens between last inspection time
            ## and lifetime set illness status to zero if the
            ## subject is right censored after last inspection time 
            ## and otherwise to user option illness.known.at.death
            if (dat$illstatus[i]){
                hit <- prodlim::sindex(eval.times=dat[i,"illtime"],jump.times=itimes,strict=TRUE)
                if (hit == (length(itimes)-1))
                    if (dat$censtime[i]<=itimes[length(itimes)])
                        out <- c(ctime,c(0,itimes)[c(1+hit,2+hit)],0)
                    else
                        out <- c(ctime,c(0,itimes)[c(1+hit,2+hit)],illness.known.at.death)
                else
                    out <- c(ctime,c(0,itimes)[c(1+hit,2+hit)],1)
                ## if (out[1]>out[2]) browser()
                out
            }
            else{## when never ill, set both interval borders to lifetime  
                c(ctime,rep(max(itimes),2),0)
            }
        }))
        colnames(interval) <- c("censtime","L","R","seen.ill")
        dat <- dat[,-c(ipos,match("censtime",names(dat)))]
        dat <- cbind(dat,interval)
        if (latent==FALSE)
            dat <- dat[,-grep("latent\\.",names(dat))]
        if (keep.inspectiontimes)
            dat$inspectiontimes <- iframe
    }
    dat$seen.exit <- 1*(dat$lifetime<dat$censtime)
    dat$observed.lifetime <- pmin(dat$lifetime,dat$censtime)
    dat$observed.illtime <- pmin(dat$illtime,dat$censtime)
    dat$L <- pmin(dat$L,dat$censtime)
    dat$R <- pmin(dat$R,dat$censtime)
    dat$observed.illtime[dat$illstatus==0] <- -9
    dat$illtime[dat$illstatus==0] <- -9
    dat
}
