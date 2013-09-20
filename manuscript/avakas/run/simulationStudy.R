library(foreach)
library(lava)
library(SmoothHazard)
library(mvtnorm)
nclust <- 200
library(doSNOW)
library("Rmpi") 
cl <- makeCluster(nclust, type = "MPI") 
registerDoSNOW(cl)
## scenario <- expand.grid(N=c(50,250,500),b01=c(0,log(2)),b02=c(0,log(2)),K=5,schedule=c(5,35),punctuality=1/20)
scenario <- expand.grid(schedule=c(0,5,20,35),N=c(50,250,500),b01=log(2),b02=log(2),K=10,punctuality=1)
runOne <- function(i,NS,seed){
    message(i)
    p <- scenario[i,,drop=FALSE]
    set.seed(seed)
    iseed <- sample(1:1000000,size=NS,replace=FALSE)
    source("~/research/SoftWare/SmoothHazard/manuscript/avakas/R/tictoc.R")
    ## setwd("~/research/SoftWare/SmoothHazard/manuscript/avakas/")
    ## source("R/tictoc.R")
    tic()
    inner <- foreach(s = 1:NS,.packages=c("SmoothHazard","lava"),.export=c("p","iseed","tic","toc")) %dopar% {
        set.seed(iseed[s])
        mod <- idmModel(K=p$K,schedule=p$schedule,punctuality=p$punctuality)
        regression(mod,from="X",to="lifetime") <- p$b02
        regression(mod,from="X",to="waittime") <- p$b02
        regression(mod,from="X",to="illtime") <- p$b01
        dat <- sim(mod,p$N)
        ## p.event <- table(dat$event)/NROW(dat)
        ## p.ill <- table(dat$ill)/NROW(dat)
        ## len.intervals <- median(dat[dat$ill==1,"R"]-dat[dat$ill==1,"L"])
        if (p$schedule==0)
            form.ill <- Hist(time=illtime,event=ill)~X
        else
            form.ill <- Hist(time=list(L,R),event=ill)~X
        ## weibull model
        try.weib <- try(weib <- idm(formula02=Hist(time=lifetime,event=status)~X,
                                    formula01=form.ill,
                                    data=dat,
                                    intensities="Weib"),silent=TRUE)
        if (inherits(try.weib,"try-error")==TRUE) {
            innerout <- list(weib=list(coef=NA,conv=NA))
        }
        else{
            innerout <- list(weib=list(coef=weib$coef,modelpar=weib$modelPar,conv=weib$converged))
        }
        ## Spline model with quantile knots 
        try.splines.quantiles <- try(suppressWarnings(splines.quantiles <- idm(formula02=Hist(time=lifetime,status)~X,
                                                                             formula01=form.ill,
                                                                             data=dat,
                                                                             maxiter=2000,
                                                                             knots="quantiles")),silent=TRUE)
        if (inherits(try.splines.quantiles,"try-error")==TRUE) {
            innerout <- c(innerout,list(splines.quantiles=list(coef=NA,conv=NA,maxiter=NA)))
        }
        else{
            innerout <- c(innerout,list(splines.quantiles=list(coef=splines.quantiles$coef,conv=splines.quantiles$converged,maxiter=splines.quantiles$niter)))
        }
        ## Spline model with equidistant knots 
        try.splines.equi <- try(suppressWarnings(splines.equi <- idm(formula02=Hist(time=lifetime,status)~X,
                                                                   formula01=form.ill,
                                                                   data=dat,
                                                                   maxiter=2000,
                                                                   knots="equidistant",
                                                                   intensities="Splines")),silent=TRUE)
        if (inherits(try.splines.equi,"try-error")==TRUE) {
            innerout <- c(innerout,list(splines.equi=list(coef=NA,conv=NA,maxiter=NA)))
        }
        else{
            innerout <- c(innerout,list(splines.equi=list(coef=splines.equi$coef,conv=splines.equi$converged,maxiter=splines.equi$niter)))
        }
        if (s %in% c(1,10,50,seq(100,1000,100))) message(paste("Simulations: ",s,", Minutes: ",round(toc()/60,2),sep=""))
        innerout
    }
    inner
}
b <- runOne(i=1,NS=1,seed=1)


## run simulation
seed <- 1735
NS <- 1000
setwd("~/research/SoftWare/SmoothHazard/manuscript/avakas/results/")
nix <- lapply(1:NROW(scenario),function(i){
    results <- runOne(i,NS=NS,seed=seed)
    save(results,file=paste("results-",i,"-seed",seed,"-NS-",NS,"-",sub(" ","-",Sys.time()),".rda",sep=""))
})
