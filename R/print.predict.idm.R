### print.predict.idm.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 16 2016 (09:50) 
## Version: 
## last-updated: Feb 27 2016 (08:41) 
##           By: Thomas Alexander Gerds
##     Update #: 30
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @export

print.predict.idm <- function(x,digits=2,...){
    lifeExpect <- is.infinite(x$t)
    cat("Predictions of an irreversible illness-death model with states (0,1,2).\n\n")
    cat("For covariate values:\n\n")
    print(x$newdata,row.names=FALSE)
    cat("\n")
    fmt <- paste0("%1.", digits[[1]], "f")
    ## nix <- lapply(names(x$transprob),function(ntp){
    ## tp <- x$transprob[[ntp]]
    ## if (length(tp)==3)
    ## paste(sprintf(tp[[1]],fmt=fmt)," [",sprintf(tp[[2]],fmt=fmt),",",sprintf(tp[[3]],fmt=fmt),"]",sep="")
    ## else
    ## sprintf(tp[[1]],fmt=fmt)
    ## })
    px <- x$transprob
    rownames(px) <- NULL
    if (lifeExpect==TRUE){
        cat("Remaining life expected sojourn times (starting at time ",x$s,"):\n\n",sep="")
        print(cbind("State at time s"=c("0","0","1"),"Expected years in states 0,1"=c("Total","In state 0","Total"),px[px$Parameter %in% c("LE.0","LE.nondiseased","LE.diseased"),]),row.names=FALSE)
    }else{
        cat("For a subject in state '0' at time ",x$s,",\npredicted state occupation probability at time ",x$t,":\n\n",sep="")
        print(cbind("State"=c(0,1,2),px[px$Parameter %in% c("p00","p01","p02"),]),row.names=FALSE)
        cat("\nThe probability p02 can be further decomposed into\ndirect and indirect transition probabilities:\n\n")
        print(cbind("Path"=c("direct","via state 1","total"),px[px$Parameter %in% c("p02_0","p02_1","p02"),]),row.names=FALSE)
        cat("\nFor a subject in state '0' at time ",x$s,",\npredicted probability of exit from state 0 until time ",x$t,":\n\n",sep="")
        print(cbind("Path"=c("via state 1","any"),px[px$Parameter %in% c("F01","F0."),]),row.names=FALSE)
        cat("\nFor a subject in state '1' at time ",x$s,",\npredicted state occupation probability at time ",x$t,":\n\n",sep="")
        print(cbind("State"=c(1,2),px[px$Parameter %in% c("p11","p12"),]),row.names=FALSE)
    }
    invisible(px)
}



#----------------------------------------------------------------------
### print.predict.idm.R ends here
