#' Paquid data set
#' 
#' Paquid data set composed of 1000 subjects selected randomly from the Paquid
#' data set of 3675 subjects.
#' 
#' 
#' @name Paq1000
#' @docType data
#' @format A data frame with 1000 rows and the following 8 columns.  \describe{
#' \item{list("dementia")}{dementia status, 0=non-demented, 1=demented}
#' \item{list("death")}{death status, 0=alive, 1=dead} \item{list("e")}{age at
#' entry in the study} \item{list("l")}{for demented subjects: age at the visit
#' before the diagnostic visit; for non-demented subjects: age at the last
#' visit (censoring age)} \item{list("r")}{for demented subjects: age at the
#' diagnostic visit; for non-demented subjects: age at the last visit
#' (censoring age)} \item{list("t")}{for dead subjects: age at death; for alive
#' subject: age at the latest news} \item{list("certif")}{primary school
#' certificate:\code{0=with certificate}, \code{1=without certificate}}
#' \item{list("gender")}{gender: \code{0=female}, \code{1=male}} }
#' @keywords datasets
#' @examples
#' 
#' data(Paq1000)
#' 
NULL





#' Plot method for an illness-death model.
#' 
#' Plot estimated baseline transition intensities from an object of class
#' \code{idmWeib} or \code{idmPl}. Pointwise confidence intervals are
#' available.
#' 
#' 
#' @param x a \code{idmWeib} class object or a \code{idmSplines} class object
#' (output from calling \code{\link{idm}} function).
#' @param conf.int logical value. Determines whether pointwise confidence
#' intervals will be plotted. The default is FALSE.
#' @param citype Type of pointwise confidence intervals, either "shadow" or
#' "lines".
#' @param add If TRUE no new plot is generated.
#' @param axes If FALSE no axes are drawn.
#' @param col Colour of the lines
#' @param lwd Thickness of the lines
#' @param lty lty of the lines
#' @param xlim Limits of the x-axis
#' @param ylim vector containing the y limites (min, max) of the plot. The
#' default is the min and the max of all the values of the three transition
#' intensities
#' @param xlab Label for the x-axis
#' @param ylab Label for the y-axis
#' @param legend Whether to show the legend.
#' @param ylim vector containing the y limites (min, max) of the plot. The
#' default is the min and the max of all the values of the three transition
#' intensities
#' @param pos.legend The location of the legend can be specified by setting
#' this argument to a single keyword from the list "bottomright", "bottom",
#' "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
#' The default is "topleft"
#' @param main title of plot
#' @param transition a vector of the transitions to plotted on separate graphs.
#' Possible value are "01","02" and "12". The default is a plot of the 3
#' intensities on the same graph.
#' @param \dots other graphical parameters.
#' @return Print a plot of the baseline transition intensities of an
#' illness-death model.
#' @author Thomas Alexander. Gerds
#' @seealso \code{\link{idm}}
#' @keywords methods
#' @examples
#' 
#' \dontrun{
#' data(Paq1000)
#' fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#' 		formula01=Hist(time=list(l,r),event=dementia)~certif,
#' 		data=Paq1000) 
#' 
#' # no pointwise confidence intervals
#' plot(fit.weib)
#' 
#' # pointwise confidence intervals for transition 0 --> 1
#' plot(fit.weib,conf.int=TRUE,transition="01")
#' }
#' 
NULL





#' Plot method for an illness-death model using a penalized likelihood
#' approach.
#' 
#' Plot estimated baseline transition intensities from a \code{idmSplines}
#' class object. Pointwise confidence intervals are available.
#' 
#' 
#' @param x a \code{idmSplines} class object (output from calling \code{idm}
#' with the option \code{intensities}="Splines".
#' @param \dots other graphical parameters like those in \code{\link{plot.idm}}
#' @return Print a plot of the baseline transition intensities of an
#' illness-death model estimated using a penalized log-likelihood approach.
#' @seealso
#' \code{\link{print.idmSplines}},\code{\link{summary.idmSplines}},\code{\link{idm}},
#' @keywords methods
#' @examples
#' 
#' 
#' 
#' \dontrun{
#' data(Paq1000)
#' 
#' fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#' 		formula01=Hist(time=list(l,r),event=dementia)~certif,
#'                 formula12=~1,
#'                 method="Splines",
#' 		data=Paq1000)
#' plot(fit.splines)
#'  
#' } 
#' 
#' 
NULL





#' Data set for survival models: right-censored and interval-censored data.
#' 
#' A simulated data frame for survival models composed of right-censored and
#' interval-censored data.
#' 
#' 
#' @name testdata
#' @docType data
#' @format A data frame with 936 observations on the following 4 variables.
#' \describe{ \item{list("l")}{for diseased subjects: left endpoint of
#' censoring interval; for non-diseased subjects: right censoring time}
#' \item{list("r")}{for diseased subjects: right endpoint of censoring
#' interval; for non-diseased subjects: right censoring time for the disease
#' event} \item{list("id")}{disease status} \item{list("cov")}{covariate} }
#' @keywords datasets
#' @examples
#' 
#' data(testdata)
#' head(testdata)
#' 
NULL



