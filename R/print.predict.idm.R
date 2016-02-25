### print.predict.idm.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 16 2016 (09:50) 
## Version: 
## last-updated: Feb 25 2016 (08:26) 
##           By: Thomas Alexander Gerds
##     Update #: 15
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @export
print.predict.idm <- function(x,digits=3,...){
    print(x$transprob,digits=digits,...)
    invisible(x$transprob)
}


#----------------------------------------------------------------------
### print.predict.idm.R ends here
