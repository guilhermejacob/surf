#' @importFrom stats model.frame model.matrix na.pass printCoefmat terms terms.formula update weights formula as.formula var pchisq chisq.test pf
#' @importFrom Matrix solve
#' @importFrom abind abind
#' @importFrom stringr str_pad
#' @importFrom scales number

#' @importFrom survey postStratify degf
#' @export
postStratify.surflow.design <- function (design, strata, population, partial = FALSE, ...) {

  {
    if (inherits(strata, "formula")) {
      mf <- substitute(model.frame(strata, data = design$variables[[1]],
                                   na.action = na.fail))
      strata <- eval.parent(mf)
    }
    strata <- as.data.frame(strata)
    sampletable <- stats::xtabs(I(1/design$prob) ~ ., data = strata)
    sampletable <- as.data.frame(sampletable)
    if (inherits(population, "table"))
      population <- as.data.frame(population)
    else if (is.data.frame(population))
      population$Freq <- as.vector(population$Freq)
    else stop("population must be a table or dataframe")
    if (!all(names(strata) %in% names(population)))
      stop("Stratifying variables don't match")
    nn <- names(population) %in% names(strata)
    if (sum(!nn) != 1)
      stop("stratifying variables don't match")
    names(population)[which(!nn)] <- "Pop.Freq"
    both <- merge(sampletable, population, by = names(strata),
                  all = TRUE)
    samplezero <- both$Freq %in% c(0, NA)
    popzero <- both$Pop.Freq %in% c(0, NA)
    both <- both[!(samplezero & popzero), ]
    if (any(onlysample <- popzero & !samplezero)) {
      print(both[onlysample, ])
      stop("Strata in sample absent from population. This Can't Happen")
    }
    if (any(onlypop <- samplezero & !popzero)) {
      if (partial) {
        both <- both[!onlypop, ]
        warning("Some strata absent from sample: ignored")
      }
      else {
        print(both[onlypop, ])
        stop("Some strata absent from sample: use partial=TRUE to ignore them.")
      }
    }
    reweight <- both$Pop.Freq/both$Freq
    both$label <- do.call("interaction", list(both[, names(strata)]))
    designlabel <- do.call("interaction", strata)
    index <- match(designlabel, both$label)
    attr(index, "oldweights") <- 1/design$prob
    design$prob <- design$prob/reweight[index]
    attr(index, "weights") <- 1/design$prob
    design$postStrata <- c(design$postStrata, list(index))
    design$call <- sys.call(-1)
    design
  }

}


# adapted functions from survey and convey packages
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
