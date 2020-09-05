#' @importFrom stats model.frame model.matrix na.pass printCoefmat terms terms.formula update weights formula as.formula var pchisq chisq.test pf confint coef qt vcov
#' @importFrom Matrix solve
#' @importFrom abind abind
#' @importFrom stringr str_pad
#' @importFrom scales number

# adapted functions from survey and convey packages
format.perc <- function (probs, digits) {
  paste(format(100 * probs, trim = TRUE, scientific = FALSE,
               digits = digits), "%")
}

h_fun <- function(incvar, w) {
  incvar <- incvar[ w >0 ]
  w <- w[ w >0 ]
  N <- sum(w)
  sd_inc <- sqrt((sum(w * incvar * incvar) - sum(w * incvar) * sum(w * incvar)/N)/N)
  h <- sd_inc/exp(0.2 * log(sum(w)))
  h
}

densfun <- function(incvar, w, x, h = NULL, FUN = "F" , na.rm=FALSE, ...) {

  if( !( FUN %in% c( "F" , "big_s" ) ) ) stop( "valid choices for `FUN=` are 'F' and 'big_s'" )
  if(na.rm){
    nas<-is.na(incvar)
    w <- w[!nas]
    incvar <- incvar[!nas]
  }
  N <- sum(w)
  if(is.null(h)) h <- h_fun(incvar,w)
  u <- (x - incvar)/h
  vectf <- exp(-(u^2)/2)/sqrt(2 * pi)
  if (FUN == "F")
    res <- sum(vectf * w)/(N * h) else {
      v <- w * incvar
      res <- sum(vectf * v)/h
    }
  res
}

ComputeQuantiles <- function( incvar , w, p , na.rm ) {

  if(na.rm){
    nas<-is.na(incvar)
    w <- w[!nas]
    incvar <- incvar[!nas]
  }

  if (any(is.na(incvar))) return(NA * p)

  if( sum( w ) == 0 ) return( NA )

  oo <- order(incvar)
  cum.w <- cumsum(w[oo])/sum(w)
  cdf <- stats::approxfun(cum.w, incvar[oo], method = "constant", f = 1, yleft = min(incvar), yright = max(incvar), ties = min)
  cdf(p)
}
