library(SmoothHazard)
library(lava)
set.seed(17)
u <- idmModel(K=10,schedule=15,punctuality=1/20)
distribution(u,"lifetime") <- coxWeibull.lvm(scale=1/100)
regression(u,from="X",to="lifetime") <- log(2)
regression(u,from="X",to="waittime") <- log(2)
testdat1 <- sim(u,50)
tmp <- idm(formula02=Hist(time=lifetime,event=event!=0)~X,
           formula01=Hist(time=list(L,R),event=ill)~X,
           data=testdat1,
           intensities="Weib")


data(Paq1000)
Paq1000$X <- rnorm(NROW(Paq1000))
fit.weib <- idm(formula02=Hist(time=t,event=death,entry=t0)~X,
                formula01=Hist(time=list(l,r),event=dementia)~X,
                data=Paq1000)
