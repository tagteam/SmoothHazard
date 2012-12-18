.First.lib <- function(lib,pkg) {
  library.dynam("SmoothHazard",pkg,lib)
}
.Last.lib <- function(lib)
  library.dynam.unload("SmoothHazard",lib)
