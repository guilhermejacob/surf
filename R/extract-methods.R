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

#' @importFrom stats confint
#' @exportS3Method  confint svymstat
confint.svymstat <- function(object , parm , level = 0.95, df = Inf, ...) {
  cf <- coef(object)
  pnames <- names(cf)
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- format.perc(a, 3)
  fac <- qt(a, df = df)
  ses <- SE(object)
  ci <- sweep( ses %o% fac , 1:2 , cf , "+" )
  dimnames( ci )[[3]] <- pct
  ci <- list( ci[,,1] , ci[,,2] )
  names( ci ) <- pct
  ci
}
