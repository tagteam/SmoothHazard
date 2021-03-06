### test-idm.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 22 2015 (13:57) 
## Version: 
## last-updated: Feb 25 2016 (08:43) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(SmoothHazard)
library(testthat)

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
 


model2 <- idm(formula02=Hist(time=t,event=death, entry=e)~1,
              formula01=Hist(time=list(l,r),event=dementia)~1,
              data=Paq1000,  method="Splines")
s=65
tseq <- 65:100
p01_sp=matrix(0,36,4)
p02_sp=matrix(0,36,4)
p00_sp=matrix(0,36,4)
for (t in tseq){
    print(t)
    p_spo <- predict(model2,s,t) 
    p01_sp[(t-s)+1,]<- c(0,p_spo$transprob[2,2],p_spo$transprob[2,3],p_spo$transprob[2,4])
    p02_sp[(t-s)+1,]<- c(0,p_spo$transprob[7,2],p_spo$transprob[7,3],p_spo$transprob[7,4])
    p00_sp[(t-s)+1,]<- c(0,p_spo$transprob[1,2],p_spo$transprob[1,3],p_spo$transprob[1,4])
}
#----------------------------------------------------------------------
### test-idm.R ends here
