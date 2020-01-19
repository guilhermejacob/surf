#' @exportS3Method coef flowstat
coef.flowstat <- function( object , ... ) object[,]

#' @exportS3Method vcov flowstat
vcov.flowstat <- function( object , ... ) attr( object , "var" )

#' @importFrom survey SE
#' @exportS3Method SE flowstat
SE.flowstat <- function( object , ... ) sqrt( attr( object , "var" ) )

#' @importFrom survey cv
#' @exportS3Method  cv flowstat
cv.flowstat <- function( object , ... ) SE.flowstat( object ) / coef.flowstat(object)
