library(SmoothHazard)
set.seed(17)
u <- idmModel(K=10,schedule=15,punctuality=1/20)
distribution(u,"lifetime") <- coxWeibull.lvm(scale=1/120000,shape=2.5)
distribution(u,"waittime") <- coxWeibull.lvm(scale=1/120000,shape=2.5)
regression(u,from="X",to="lifetime") <- log(2)
regression(u,from="X",to="waittime") <- log(2)
testdat1 <- sim(u,50)
## does not converge
system.time(stmp1 <- idm(formula02=Hist(time=lifetime,event=status)~X,
                         formula01=Hist(time=illtime,event=ill)~X,
                         data=testdat1,
                         n.knots=c(7,7,7),
                         intensities="Splines"))
## gives a segmentation fault when run twice in a row
system.time(stmp1 <- idm(formula02=Hist(time=lifetime,event=status)~X,
                         formula01=Hist(time=illtime,event=ill)~X,
                         data=testdat1,
                         n.knots=c(2,2,2),
                         intensities="Splines"))
system.time(stmp1 <- idm(formula02=Hist(time=lifetime,event=status)~X,
                         formula01=Hist(time=illtime,event=ill)~X,
                         data=testdat1,
                         n.knots=c(2,2,2),
                         intensities="Splines"))
