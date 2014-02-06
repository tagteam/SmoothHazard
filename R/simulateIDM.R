##' Function to generate an illness-death model for simulation.
##'
##' Based on the functionality of the lava PACKAGE a latent variable model is
##' generated which contains the variables illtime, lifetime, waiting time and
##' censoring time.
##' @title Generate illness-death model objects
##' @param scale Rate for Weibull models 
##' @param cens Type of censoring 
##' @param K Number of intervals
##' @param schedule Average time for the next inspection time.
##' @param punctuality Standard deviation of time for the next inspection time.
##' @return A latent variable model object \code{lvm}
##' @author Thomas Alexander Gerds
idmModel <- function(scale=1/100,
                     cens="interval",
                     K=5,
                     schedule=1/scale,
                     punctuality=10*scale){
    require(lava)
    ## ===============================================
    ## illness-death-model
    ## ===============================================
    ## model waiting time in state 0
    ## based on latent "ill" and "life" times
    ## separately model waiting time in state "ill"
    ## ===============================================
    idm <- lvm()
    latent(idm) <- c("illtime","lifetime","waittime")
    distribution(idm,"illtime") <- coxWeibull.lvm(scale=scale)
    distribution(idm,"lifetime") <- coxWeibull.lvm(scale=scale)
    idm <- eventTime(idm,time~min(illtime=1,lifetime=2),"event")
    distribution(idm,"waittime") <- coxWeibull.lvm(scale=scale)
    ## ===============================================
    ## censoring model
    ## ===============================================
    ## if (cens=="right"){
    distribution(idm,"censtime") <- coxWeibull.lvm(scale=scale)
    ## } else{
    if (cens=="interval"){
        for (k in 1:K){
            distribution(idm,paste("inspection",k,sep="")) <- normal.lvm(mean=schedule,sd=punctuality)
        }
    }
    ## }
    class(idm) <- c("idmModel","lvm")
    attr(idm,"cens") <- cens
    idm
}
##' Function to simulate illness-death model data
##'
##' Based on the functionality of the lava PACKAGE 
##' @title Simulate illness-death model data
##' @param x An \code{idmModel} object as obtained with
##' \code{idmModel}
##' @param n Number of observations
##' @param compliance Probability of missing inspection time.
##' @param ... Extra arguments given to \code{sim}
##' @return A data set with interval censored observations from an illness-death model
##' @author Thomas Alexander Gerds
sim.idmModel <- function(x,
                         n=100,
                         compliance=1,
                         ...){
    # simulate latent data
    class(x) <- "lvm"
    dat <- sim(x,n=n,...)
    # construct lifetime for ill subjects
    ill <- dat$event==1
    dat$lifetime[ill] <- dat$illtime[ill]+dat$waittime[ill]
    # reset illtime for subjects that were never ill 
    dat$illtime[!ill] <- dat$lifetime[!ill]
    cens <- attr(x,"cens")
    if (cens=="interval") {
        ipos <- grep("inspection[0-9]+",names(dat))
        interval <- do.call("rbind",lapply(1:n,function(i){
            ## cumulate times between inspections
            itimes <- unique(cumsum(c(0,pmax(0,dat[i,ipos,drop=TRUE]))))
            itimes <- c(itimes[itimes<dat$lifetime[i]],dat[i,"lifetime"])
            if (compliance!=1){
                comp <- rbinom(length(itimes),1,compliance)
                unique(round(itimes[comp==1]))
            }
            if (ill[i]){
                hit <- sindex(eval.times=dat[i,"illtime"],jump.times=itimes,strict=TRUE)
                c(c(0,itimes)[c(1+hit,2+hit)],1)
            }
            else{
                c(rep(max(itimes),2),0)
            }
        }))
        colnames(interval) <- c("L","R","ill")
        dat <- cbind(dat[,-c(ipos,match(c("time","waittime"),names(dat)))],interval)
        ## right censored?
        dat$censtime <- pmax(dat$censtime,dat$R)
    }
    dat$status <- 1*(dat$lifetime<dat$censtime) 
    dat$lifetime <- pmin(dat$lifetime,dat$censtime)
    ## dat <- dat[,-match(c("censtime","illtime"),names(dat))]
    dat
}

