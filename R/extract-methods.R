#' @exportS3Method coef flowstat
coef.flowstat <- function( object , ... ) {
  attr(object, "statistic") <- NULL
  attr(object, "deff") <- NULL
  attr(object, "var") <- NULL
  attr(object, "rounds") <- NULL
  attr(object, "formula") <- NULL
  unclass(object)
}

#' @exportS3Method vcov flowstat
vcov.flowstat <- function( object , ... ) unclass( attr( object , "var" ) )

#' @importFrom survey SE
#' @exportS3Method SE flowstat
SE.flowstat <- function( object , ... ) unclass( sqrt( attr( object , "var" ) ) )

#' @importFrom survey cv
#' @exportS3Method  cv flowstat
cv.flowstat <- function( object , ... ) SE.flowstat( object ) / coef.flowstat(object)
