#' Subset of surflow
#'
#' Restrict a survey flow design to a subpopulation, keeping the original survey flow design information about number of clusters, strata.
#'
#' @param x  a surflow design object
#' @param subset	 An expression specifying the subpopulation
#' @param rounds  a vector of integers indicating which round to apply subsetting condition.
#' Defaults to \code{rounds = NULL}, and the expression is applied to all rounds, keeping observation that were true in at least of the rounds.
#' @param on.full.design  A logical, indicating wether the subset should be applied to full design weights. Only meaningful for endogenous categories.
#' @param ...  future expansion
#'
#' @examples
#' # load data
#' data( "real" )
#'
#' # create surflow design object
#' flowdes <-
#'   sfydesign( ids = ~ upa ,
#'              probs = ~ v1027 ,
#'              strata = ~ estrato ,
#'              data = list( pnadc1 , pnadc2 ) ,
#'              nest = TRUE )
#'
#' # subset
#' flowdes <- subset( flowdes , v2009 >= 14 )
#'
#' @method subset surflow.design
#' @export
subset.surflow.design <-
  function (x , subset , rounds = NULL , ... , on.full.design = FALSE ) {

    if ( is.null( rounds ) ) rounds <- seq_along( x$variables )
    if ( is.numeric( rounds ) ) {
      rounds <- unique( as.integer( rounds ) )
      if ( !all( rounds %in% seq_along( x$variables ) ) ) stop( "invalid index; check ?subset.surflow.design for examples.")
    } else stop( "invalid index; check ?subset.surflow.design for examples.")

    # apply subsetting
    e <- substitute( subset )
    r <- rep( TRUE , nrow( x$variables[[1]] ) )
    for ( z in rounds ) {
      s <- eval(e, x$variables[[z]], parent.frame())
      # r <- r & s & !is.na(s) # original line
      r <- r & ( s | is.na(s) )
    }

    if ( "svyrep.design" %in% class(x) ) {
      if ( on.full.design ) {
        x$fullpweights[!r] <- 0
        x$fullrepweights$weights[!r,] <- 0
      }
      x$pweights[!r] <- 0
      x$repweights$weights[!r,] <- 0
      x$degf <- NULL
      x$degf <- survey::degf(x)
    } else {
      if ( on.full.design ) x$fullprob[!r] <- Inf
      x$prob[!r] <- Inf
    }
    return( x )
  }

