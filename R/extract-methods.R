#' @exportS3Method coef svymstat
coef.svymstat <- function( object , ... ) {
  attr(object, "categories") -> these.dimnames
  object <- unclass(object)
  attr(object, "statistic") <- NULL
  attr(object, "var") <- NULL
  attr(object, "categories") <- NULL
  object <- matrix( object , nrow = sapply( these.dimnames , length )[1] , dimnames = these.dimnames , byrow = TRUE )
  object
}

#' @exportS3Method vcov svymstat
vcov.svymstat <- function( object , ... ) {
  var.mat <- attr( object , "var" )
  var.mat
}

#' @importFrom survey SE
#' @exportS3Method SE svymstat
SE.svymstat <- function( object , ... ) {
  var.mat <- attr( object , "var" )
  var.mat <- matrix( diag( var.mat ) , nrow = length( attr( object , "categories" )[[1]] ) , byrow = TRUE )
  dimnames( var.mat ) <- attr( object , "categories" )
  sqrt( var.mat )
}

#' @importFrom survey cv
#' @exportS3Method  cv svymstat
cv.svymstat <- function( object , ... ) {
  se.mat   <- SE( object )
  coef.mat <- coef( object )
  cv.mat <- se.mat / coef.mat
  cv.mat
}

#' @importFrom stats confint
#' @exportS3Method  confint svymstat
confint.svymstat <- function(object , parm , level = 0.95, df = Inf, ...) {
  coef.mat <- invisible( coef( object ) )
  pnames <- names( coef.mat )
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- format.perc(a, 3)
  fac <- qt(a, df = df)
  ses <- SE( object )
  ci.mat <- sweep( ses %o% fac , 1:2 , coef.mat , "+" )
  dimnames( ci.mat )[[3]] <- pct
  ci.mat <- list( ci.mat[,,1] , ci.mat[,,2] )
  names( ci.mat ) <- pct
  ci.mat
}

#' @importFrom survey svycontrast
#' @export
svycontrast.svymstat <- function( stat , contrasts , ... ) {

  # remove the "svymstat" class from the current object
  class( stat ) <- setdiff( class(stat) , "svymstat" )

  # apply usual function
  survey::svycontrast( stat , contrasts , ... )

}
