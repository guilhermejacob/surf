#' Subset of surflow
#'
#' Restrict a survey flow design to a subpopulation, keeping the original survey flow design information about number of clusters, strata.
#'
#' @param x  a surflow design object
#' @param subset	 An expression specifying the subpopulation
#' @param rounds  a vector of integers indicating which round to apply subsetting condition. Defaults to \code{rounds = 0}, and the expression is applied to the first round.
#' @param on.full.design  A logical, indicating wether the subset should be applied to full design weights. Only meaningful for endogenous categories.
#' @param ...  future expansion
#'
#' @examples
#' # load data
#' data( "artificial" )
#'
#' # create surflow design object
#' flowdes <-
#'   sfydesign( ids = ~ upa ,
#'              probs = ~ longprob ,
#'              strata = ~ estrato ,
#'              data = list( dfr0 , dfr1 ) ,
#'              nest = TRUE )
#'
#' # subset
#' flowdes <- subset( flowdes , v2007 == 1 )
#'
#' @method subset surflow.design
#' @export
subset.surflow.design <-
  function (x , subset , rounds = 0 , on.full.design = FALSE , ... ) {

    if ( is.null( rounds ) ) rounds <- seq_along( x$variables ) - 1
    if ( is.numeric( rounds ) ) {
      rounds <- unique( as.integer( rounds) )
      if ( !all( rounds %in% c( seq_along( x$variables ) - 1 ) ) ) stop( "invalid index; check ?subset.surflow.design for examples.")
    } else stop( "invalid index; check ?subset.surflow.design for examples.")

    if ( "svyrep.design" %in% class(x) ) {
      e <- substitute( subset )
      r <- rep( TRUE , nrow( x$variables[[1]] ) )
      for ( z in rounds ) {
        s <- eval(e, x$variables[[z+1]], parent.frame())
        # r <- r & s & !is.na(s) # original line
        r <- r & ( s | is.na(s) )
      }
      if ( on.full.design ) {
        x$fullpweights[!r] <- 0
        x$fullrepweights$weights[!r,] <- 0
      }
      x$pweights[!r] <- 0
      x$repweights$weights[!r,] <- 0
      x$degf <- NULL
      x$degf <- survey::degf(x)
    } else {
      e <- substitute( subset )
      r <- rep( TRUE , nrow( x$variables[[1]] ) )
      for ( z in rounds ) {
        s <- eval(e, x$variables[[z+1]], parent.frame() )
        # r <- r & s & !is.na(s) # original line
        r <- r & ( s | is.na(s) )
      }
      if ( on.full.design ) x$fullprob[!r] <- Inf
      x$prob[!r] <- Inf
    }
    return( x )

  }

