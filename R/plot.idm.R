plot.idmWeib <- plot.idmPl <- plot.idm <- function(x,conf.int=FALSE, 
#ylim=FALSE,type.plot="intensity", 
ylim,pos.legend="topleft", main, transition, ...){ 

	if(missing(main)) main<-""
	if(!missing(transition)){
		if(sum(c("01","02","12") %in% transition)!=length(transition))stop("transition argument must be in 'c('01','02','12')'")
	}

	if(class(x) == "idmWeib"){
	
			if (missing(transition)){
				vec <- c("transition intensity function for 0--->1", 
					"transition intensity function for 0--->2",
					"transition intensity function for 1--->2")

                                if (missing(ylim)){
					ylim <- c(min(x$hazard01,x$hazard02,x$hazard12),max(x$hazard01,x$hazard02,x$hazard12))
				}
					
				if(conf.int){
					matplot(x$time, cbind(x$hazard01,x$lowerHazard01,x$upperHazard01),type="l",lty=c(1,1,1),col=c(1,1,1),lwd=c(2,1,1), xlab="Time",ylab="Weibull transition intensities",ylim=ylim, main=main)
					matlines(x$time, cbind(x$hazard02,x$lowerHazard02,x$upperHazard02),ylim=ylim, type="l",lwd=c(2,1,1),lty=c(2,2,2),col=c(2,2,2))
					matlines(x$time, cbind(x$hazard12,x$lowerHazard12,x$upperHazard12),ylim=ylim, type="l",lwd=c(2,1,1),lty=c(3,3,3),col=c(3,3,3))
				}else{
					matplot(x$time, cbind(x$hazard01,x$hazard02,x$hazard12),type="l",lwd=2,lty=c(1,2,3),col=c(1,2,3), xlab="Time",ylab="Weibull transition intensities",ylim=ylim, main=main)
				}	
				legend(pos.legend,c("01","02","12"),lty=c(1,2,3),lwd=c(2,2,2),col=c(1,2,3))
			}else{
				if(length(transition)==3){
					par(mfrow=c(3,1))

					if("01" %in% transition){
					        if (missing(ylim)){
							ylim <- c(min(c(x$hazard01,x$lowerHazard01,x$upperHazard01)),max(c(x$hazard01,x$lowerHazard01,x$upperHazard01)))
						}
						if(!conf.int){
							plot(x$time, x$hazard01,type="l",col=1,lwd=2, xlab="Time",ylab="Weibull transition intensity 01",main=main)
						}else{
							matplot(x$time, cbind(x$hazard01,x$lowerHazard01,x$upperHazard01), col=c(1,1,1),type="l",lty=c(1,1,1),lwd=c(2,1,1),ylim=ylim
							,xlab="Time",ylab="Weibull transition intensity 01", main=main)	
#							legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
						}
						
					}
					if("02"%in% transition){
						if (missing(ylim)){
							ylim <- c(min(c(x$hazard02,x$lowerHazard02,x$upperHazard02)),max(c(x$hazard02,x$lowerHazard02,x$upperHazard02)))
						}
						if(!conf.int){
							plot(x$time, x$hazard02, type="l",col=2,lwd=2,lty=2, xlab="Time",ylab="Weibull transition intensity 02",main=main)
						}else{
							matplot(x$time, cbind(x$hazard02,x$lowerHazard02,x$upperHazard02),type="l",col=c(2,2,2),lwd=c(2,1,1),lty=c(2,2,2),ylim=ylim
							,xlab="Time",ylab="Weibull transition intensity 02",main=main)	
#							legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
						}
					}				
					if("12" %in% transition){
						if (missing(ylim)){
							ylim <- c(min(c(x$hazard12,x$lowerHazard12,x$upperHazard12)),max(c(x$hazard12,x$lowerHazard12,x$upperHazard12)))
						}
						if(!conf.int){
							plot(x$time, x$hazard12, type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Weibull transition intensity 12", main=main)
						}else{
							matplot(x$time, cbind(x$hazard12,x$lowerHazard12,x$upperHazard12),type="l",col=c(3,3,3),lwd=c(2,1,1),lty=c(3,3,3),ylim=ylim
							 ,xlab="Time",ylab="Weibull transition intensity 12",main=main)	
#							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
						}
					}
				}else{	
				if(length(transition)==2){
					par(mfrow=c(2,1))

					if("01" %in% transition){
					        if (missing(ylim)){
							ylim <- c(min(x$hazard01,x$lowerHazard01,x$upperHazard01),max(x$hazard01,x$lowerHazard01,x$upperHazard01))
						}
						if(!conf.int){
							plot(x$time, x$hazard01,type="l",col=1,lwd=2, xlab="Time",ylab="Weibull transition intensity 01",main=main)
						}else{
							matplot(x$time, cbind(x$hazard01,x$lowerHazard01,x$upperHazard01), col=c(1,1,1),type="l",lty=c(1,1,1),lwd=c(2,1,1),ylim=ylim
							,xlab="Time",ylab="Weibull transition intensity 01", main=main)	
#							legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
						}
						
					}
					if("02"%in% transition){
						if (missing(ylim)){
							ylim <- c(min(c(x$hazard02,x$lowerHazard02,x$upperHazard02)),max(c(x$hazard02,x$lowerHazard02,x$upperHazard02)))
						}
						if(!conf.int){
							plot(x$time, x$hazard02, type="l",col=2,lwd=2,lty=2, xlab="Time",ylab="Weibull transition intensity 02",main=main)
						}else{
							matplot(x$time, cbind(x$hazard02,x$lowerHazard02,x$upperHazard02),type="l",col=c(2,2,2),lwd=c(2,1,1),lty=c(2,2,2),ylim=ylim
							,xlab="Time",ylab="Weibull transition intensity 02",main=main)	
#							legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
						}
					}				
					if("12" %in% transition){
						if (missing(ylim)){
							ylim <- c(min(c(x$hazard12,x$lowerHazard12,x$upperHazard12)),max(c(x$hazard12,x$lowerHazard12,x$upperHazard12)))
						}
						if(!conf.int){
							plot(x$time, x$hazard12, type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Weibull transition intensity 12", main=main)
						}else{
							matplot(x$time, cbind(x$hazard12,x$lowerHazard12,x$upperHazard12),type="l",col=c(3,3,3),lwd=c(2,1,1),lty=c(3,3,3),ylim=ylim
							 ,xlab="Time",ylab="Weibull transition intensity 12",main=main)	
#							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
						}
					}
				}else{
					if(transition=="01"){
						if (missing(ylim)){
							ylim <- c(min(x$hazard01,x$lowerHazard01,x$upperHazard01),max(x$hazard01,x$lowerHazard01,x$upperHazard01))
						}
						if(!conf.int){
							plot(x$time,x$hazard01,type="l",col=1,lwd=2,xlab="Time",ylab="Weibull transition intensity 01",main=main)
						}else{
							matplot(x$time, cbind(x$hazard01,x$lowerHazard01,x$upperHazard01),type="l",col=c(1,1,1),lty=c(1,1,1),lwd=c(2,1,1),ylim=ylim
							 ,xlab="Time",ylab="Weibull transition intensity 01", main=main)
#							 legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,1,1),cex=0.7)
						}
					}
					if(transition=="02"){
						if (missing(ylim)){
							ylim <- c(min(x$hazard02,x$lowerHazard02,x$upperHazard02),max(x$hazard02,x$lowerHazard02,x$upperHazard02))
						}
						if(!conf.int){
							plot(x$time, x$hazard02,type="l",col=2,lwd=2,lty=2,xlab="Time",ylab="Weibull transition intensity 02",main=main)
						}else{
							matplot(x$time,cbind(x$hazard02,x$lowerHazard02,x$upperHazard02),type="l",col=c(2,2,2),lwd=c(2,1,1),lty=c(2,2,2),ylim=ylim
							 ,xlab="Time",ylab="Weibull transition intensity 02",main=main)	
#							 legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(2,2,2),cex=0.7)	
						}
					}				
					if(transition=="12"){
						if (missing(ylim)){
							ylim <- c(min(x$hazard12,x$lowerHazard12,x$upperHazard12),max(x$hazard12,x$lowerHazard12,x$upperHazard12))
						}
						if(!conf.int){
							plot(x$time, x$hazard12, type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Weibull transition intensity 12", main=main)
						}else{
							matplot(x$time, cbind(x$hazard12,x$lowerHazard12,x$upperHazard12),type="l",col=c(3,3,3),lwd=c(2,1,1),lty=c(3,3,3),ylim=ylim
							 ,xlab="Time",ylab="Weibull transition intensity 12", main=main)	
#							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(3,3,3),cex=0.7)
						}
					}
				}
			}	
		}
	}else{
		if (missing(transition)){
				vec <- c("transition intensity function for 0--->1", 
					"transition intensity function for 0--->2",
					"transition intensity function for 1--->2")
				xmin <- min(x$time)
				xmax <- max(x$time)
				if (missing(ylim)){
					ymin <- min(x$hazard01,x$hazard02,x$hazard12)
					ymax <- max(x$hazard01,x$hazard02,x$hazard12)
					ymin01 <- min(x$hazard01,x$lowerHazard01,x$upperHazard01)
					ymax01 <- max(x$hazard01,x$lowerHazard01,x$upperHazard01)
					ymin02 <- min(x$hazard02,x$lowerHazard02,x$lowerHazard02)
					ymax02 <- max(x$hazard02,x$lowerHazard02,x$lowerHazard02)
					ymin12 <- min(x$hazard12,x$lowerHazard12,x$lowerHazard12)
					ymax12 <- max(x$hazard12,x$lowerHazard12,x$lowerHazard12)
					ylim <- c(ymin,ymax)
                                }
                                
				if(conf.int){
					matplot(x$time[,1], cbind(x$hazard01,x$lowerHazard01,x$upperHazard01),type="l",lwd=c(2,1,1),lty=c(1,1,1),col=c(1,1,1),
					xlab="Time",ylab="Splines transition intensities",xlim=c(xmin,xmax),ylim=ylim, main=main)
					matlines(x$time[,2], cbind(x$hazard02,x$lowerHazard02,x$upperHazard02),xlim=c(xmin,xmax),ylim=ylim,type="l",lty=c(2,2,2),lwd=c(2,1,1),col=c(2,2,2))
					matlines(x$time[,3], cbind(x$hazard12,x$lowerHazard12,x$upperHazard12),xlim=c(xmin,xmax),ylim=ylim,type="l",lty=c(3,3,3),lwd=c(2,1,1),col=c(3,3,3))
				}else{
					plot(x$time[,1], x$hazard01, col=1, type="l",lwd=2, xlab="Time",ylab="Splines transition intensities",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
					 main=main,)
					lines(x$time[,2], x$hazard02, col=2, type="l",lty=2,lwd=2,xlim=c(xmin,xmax),ylim=ylim)
					lines(x$time[,3], x$hazard12, col=3, type="l",lty=3,lwd=2,xlim=c(xmin,xmax),ylim=ylim)
				}
				legend(pos.legend,c("01","02","12"),lty=c(1,2,3),lwd=c(2,2,2),col=c(1,2,3))	
			}else{
				if(length(transition)==3){
					par(mfrow=c(3,1))
					if("01" %in% transition){
						if (missing(ylim)){
							ylim <- c(min(x$hazard01,x$lowerHazard01,x$upperHazard01),max(x$hazard01,x$lowerHazard01,x$upperHazard01))
						}
						if(!conf.int){
							plot(x$time[,1],x$hazard01,type="l",col=1,lwd=2, xlab="Time",ylab="Splines transition intensity 01", main=main)
						}else{
							matplot(x$time[,1],cbind(x$hazard01,x$lowerHazard01,x$upperHazard01),type="l",col=c(1,1,1),lwd=c(2,1,1),lty=c(1,1,1),ylim=ylim
							 ,xlab="Time",ylab="Splines transition intensity", main=main)
#							 legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
						}
					}
					if("02"%in% transition){
						if (missing(ylim)){
							ylim <- c(min(c(x$hazard02,x$lowerHazard02,x$upperHazard02)),max(c(x$hazard02,x$lowerHazard02,x$upperHazard02)))
						}
						if(!conf.int){
							plot(x$time[,2],x$hazard02,type="l",col=2,lwd=2,lty=2, xlab="Time",ylab="Splines transition intensity 02", main=main)
						}else{
							matplot(x$time[,2],cbind(x$hazard02,x$lowerHazard02,x$upperHazard02),type="l",col=c(2,2,2),lty=c(2,2,2),lwd=c(2,1,1),ylim=ylim
							 ,xlab="Time",ylab="Splines transition intensity", main=main)
#							 legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(2,2,2),cex=0.7)
						}
					}				
					if("12" %in% transition){
						if (missing(ylim)){
							ylim <- c(min(c(x$hazard12,x$lowerHazard12,x$upperHazard12)),max(c(x$hazard12,x$lowerHazard12,x$upperHazard12)))
						}
						if(!conf.int){
							plot(x$time[,3],x$hazard12,type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Splines transition intensity 12", main=main)
						}else{
							matplot(x$time[,3],cbind(x$hazard12,x$lowerHazard12,x$upperHazard12),type="l",col=c(3,3,3),lty=c(3,3,3),lwd=c(2,1,1),ylim=ylim
							 ,xlab="Time",ylab="Splines transition intensity", main=main)
#							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(3,3,3),cex=0.7)
						}
					}
			}else{
				if(length(transition)==2){
					par(mfrow=c(2,1))
					if("01" %in% transition){
						if (missing(ylim)){
							ylim <- c(min(x$hazard01,x$lowerHazard01,x$upperHazard01),max(x$hazard01,x$lowerHazard01,x$upperHazard01))
						}
						if(!conf.int){
							plot(x$time[,1],x$hazard01, type="l",col=1,lwd=2, xlab="Time",ylab="Splines transition intensity 01", main=main)
						}else{
							matplot(x$time[,1],cbind(x$hazard01,x$lowerHazard01,x$upperHazard01),type="l",col=c(1,1,1),lwd=c(2,1,1),lty=c(1,1,1),ylim=ylim
							 ,xlab="Time",ylab="Splines transition intensity", main=main)
#							 legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,1,1),cex=0.7)
						}
					}
					if("02"%in% transition){
						if (missing(ylim)){
							ylim <- c(min(x$hazard02,x$lowerHazard02,x$upperHazard02),max(x$hazard02,x$lowerHazard02,x$upperHazard02))
						}
						if(!conf.int){
							plot(x$time[,2],x$hazard02, type="l",col=2,lwd=2,lty=2, xlab="Time",ylab="Splines transition intensity 02", main=main)
						}else{
							matplot(x$time[,2],cbind(x$hazard02,x$lowerHazard02,x$upperHazard02),type="l",col=c(2,2,2),lty=c(2,2,2),lwd=c(2,1,1),ylim=ylim
							 ,xlab="Time",ylab="Splines transition intensity", main=main)
#							 legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(2,2,2),cex=0.7)
						}
					}				
					if("12" %in% transition){
						if (missing(ylim)){
							ylim <- c(min(x$hazard12,x$lowerHazard12,x$upperHazard12),max(x$hazard12,x$lowerHazard12,x$upperHazard12))
						}
						if(!conf.int){
							plot(x$time[,3],x$hazard12, type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Splines transition intensity 12", main=main)
						}else{
							matplot(x$time[,3],cbind(x$hazard12,x$lowerHazard12,x$upperHazard12),type="l",col=c(3,3,3),lty=c(3,3,3),lwd=c(2,1,1),ylim=ylim
							 ,xlab="Time",ylab="Splines transition intensity", main=main)
#							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(3,3,3),cex=0.7)
						}
					}					
				}else{
					if(transition=="01"){
						if (missing(ylim)){
							ylim <- c(min(x$hazard01,x$lowerHazard01,x$upperHazard01),max(x$hazard01,x$lowerHazard01,x$upperHazard01))
						}
						if(!conf.int){
							plot(x$time[,1],x$hazard01,type="l",col=1,lwd=2,xlab="Time",ylab="Splines transition intensity 01",main=main)
						}else{
							matplot(x$time[,1],cbind(x$hazard01,x$lowerHazard01,x$upperHazard01),type="l",col=c(1,1,1),lwd=c(2,1,1),lty=c(1,1,1),ylim=ylim
							 ,xlab="Time",ylab="Splines transition intensity", main=main)
#							 legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,1,1),cex=0.7)
						}
					}
					if(transition=="02"){
						if (missing(ylim)){
							ylim <- c(min(x$hazard02,x$lowerHazard02,x$upperHazard02),max(x$hazard02,x$lowerHazard02,x$upperHazard02))
						}
						if(!conf.int){
							plot(x$time[,2],x$hazard02, type="l",col=2,lwd=2,lty=2, xlab="Time",ylab="Splines transition intensity 02", main=main)
						}else{
							matplot(x$time[,2],cbind(x$hazard02,x$lowerHazard02,x$upperHazard02),type="l",col=c(2,2,2),lty=c(2,2,2),lwd=c(2,1,1),ylim=ylim
							 ,xlab="Time",ylab="Splines transition intensity", main=main)
#							 legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(2,2,2),cex=0.7)
						}
					}				
					if(transition=="12"){
						if (missing(ylim)){
							ylim <- c(min(x$hazard02,x$lowerHazard02,x$upperHazard02),max(x$hazard02,x$lowerHazard02,x$upperHazard02))
						}
						if(!conf.int){
							plot(x$time[,3],x$hazard12, type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Splines transition intensity 12", main=main)
						}else{
							matplot(x$time[,3],cbind(x$hazard12,x$lowerHazard12,x$upperHazard12),type="l",col=c(3,3,3),lty=c(3,3,3),lwd=c(2,1,1),ylim=ylim
							 ,xlab="Time",ylab="Splines transition intensity", main=main)
#							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(3,3,3),cex=0.7)
						}
					}
				}
			}
		}
	}
	return(invisible())
}
