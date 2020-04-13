#' @exportS3Method coef svymstat
coef.svymstat <- function( object , ... ) {
  object <- unclass(object)
  attr(object, "statistic") <- NULL
  attr(object, "var") <- NULL
  object
}

#' @exportS3Method vcov svymstat
vcov.svymstat <- function( object , ... ) {
  object <- attr(object, "var")
  object
}

#' @importFrom survey SE
#' @exportS3Method SE svymstat
SE.svymstat <- function( object , ... ) {
  vmat <- attr(object, "var")
  sqrt( vmat )
}

#' @importFrom survey cv
#' @exportS3Method  cv svymstat
cv.svymstat <- function( object , ... ) {
  cmat <- unclass(object)
  attr(cmat, "statistic") <- NULL
  attr(cmat, "var") <- NULL
  vmat <- attr(object, "var")
  semat <- sqrt(vmat)
  semat / cmat
}
