### test-idm.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 22 2015 (13:57) 
## Version: 
## last-updated: Oct 22 2015 (16:35) 
##           By: Thomas Alexander Gerds
##     Update #: 1
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
tmp <- Paq1000
tmp$t <- 200
fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
                formula01=Hist(time=list(l,r),event=dementia)~certif,
                data=tmp)

#----------------------------------------------------------------------
### test-idm.R ends here
