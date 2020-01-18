##' @importFrom stats coef
##' @export
stats::coef
##' @method coef flowstat
##' @export
coef.flowstat <- function( object , ... ) object[,]

##' @importFrom stats vcov
##' @export
stats::vcov
##' @method vcov flowstat
##' @export
vcov.flowstat <- function( object , ... ) attr( object , "var" )

# flowstat SE method
##' @importFrom survey SE
##' @export
survey::SE
##' @method SE flowstat
##' @export
SE.flowstat <- function( object , ... ) sqrt( attr( object , "var" ) )

# flowstat cv method
##' @importFrom survey cv
##' @export
survey::cv
##' @method cv flowstat
##' @export
cv.flowstat <- function( object , ... ) SE.flowstat( object ) / coef.flowstat(object)
