### print.predict.idm.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 16 2016 (09:50) 
## Version: 
## last-updated: Feb 16 2016 (16:03) 
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
#' @export
print.predict.idm <- function(x,...){
    tab <- cbind(Parameter=names(x$transprob),data.frame(matrix(unlist(x$transprob),ncol=3,byrow=TRUE)))
    colnames(tab)[2:4] <- c("Estimate","lower.95","upper.95")
    print(tab,...)
    invisible(tab)
}


#----------------------------------------------------------------------
### print.predict.idm.R ends here
