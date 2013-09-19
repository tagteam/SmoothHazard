idmModel <- function(scale=1/100,cens="interval",K=5,schedule=1/scale,punctuality=10*scale){
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
    distribution(idm,"illtime") <- coxWeibull.lvm(scale=1/100)
    distribution(idm,"lifetime") <- coxWeibull.lvm(scale=1/100)
    idm <- eventTime(idm,time~min(illtime=1,lifetime=2),"event")
    distribution(idm,"waittime") <- coxWeibull.lvm(scale=1/100)
    ## ===============================================
    ## censoring model
    ## ===============================================
    ## if (cens=="right"){
    distribution(idm,"censtime") <- coxWeibull.lvm(scale=1/100)
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
sim.idmModel <- function(x,n=100,compliance=1,p=NULL,normal=FALSE,cond=FALSE,sigma = 1, rho=0.5,X,unlink=FALSE,...){
    # simulate latent data
    class(x) <- "lvm"
    dat <- sim(x,n=n,...)
    # construct lifetime for ill subjects
    ill <- dat$event==1
    dat$lifetime[ill] <- dat$illtime[ill]+dat$waittime[ill]
    # remove illtime for never ill subjects
    dat$illtime[!ill] <- NA
    cens <- attr(x,"cens")
    ## if (cens=="right"){
    # right censoring
    ## } else
    if (cens=="interval") {
        ipos <- grep("inspection[0-9]+",names(dat))
        interval <- do.call("rbind",lapply(1:n,function(i){
            ## cumulate times between inspections
            itimes <- cumsum(c(0,pmax(0,dat[i,ipos,drop=TRUE])))
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
    }
    ## right censored?
    dat$censtime <- pmax(dat$censtime,dat$R)
    dat$status <- 1*(dat$lifetime<dat$censtime)
    dat$lifetime <- pmin(dat$lifetime,dat$censtime)
    ## dat <- dat[,-match(c("censtime","illtime"),names(dat))]
    dat
}

