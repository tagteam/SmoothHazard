### test-idm.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 22 2015 (13:57) 
## Version: 
## last-updated: Feb 25 2016 (08:43) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(SmoothHazard)
data(Paq1000)
fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
                formula01=Hist(time=list(l,r),event=dementia)~certif,
                data=Paq1000)

library(SmoothHazard)
data(Paq1000)
fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
                formula01=Hist(time=list(l,r),event=dementia)~certif,
                formula12 = ~ 1,
                data=Paq1000)
pred <- predict(fit.weib,70,t=80,newdata=data.frame(certif=1), conf.int = TRUE, nsim=4)
pred <- predict(fit.weib,70,t=80,newdata=data.frame(certif=1), conf.int = TRUE, nsim=4,lifeExpect=TRUE)
pred
 

#----------------------------------------------------------------------
### test-idm.R ends here
