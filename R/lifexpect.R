lifexpect <- function(object,s,Z01,Z02,Z12,nsim=1000,CI=TRUE,...) {

  xx <- object
  if (inherits(xx,"idmSplines")) 
    return(lifexpect.idmSplines(xx,s,Z01,Z02,Z12,nsim,CI,...))
    else
      return(lifexpect.idmWeib(xx,s,Z01,Z02,Z12,nsim,CI,...))

}
