#' Subset of survflow
#'
#' Restrict a survey flow design to a subpopulation, keeping the original survey flow design information about number of clusters, strata.
#'
#' @method subset survflow.design
#'
#' @param x  a survflow design object
#' @param subset	 An expression specifying the subpopulation
#' @param rounds  a vector of integers indicating which round to apply subsetting condition. Defaults to \code{rounds = 0}, and the expression is applied to the first round.
#'
#' @examples
#' # load data
#' data( "initial" )
#' data( "final" )
#'
#' # create survflow design object
#' flowdes <-
#'   svyflowdesign( ids = ~ upa ,
#'                  strata = ~ estrato ,
#'                  probs = ~ longprob ,
#'                  data.list = list( f.qtr , s.qtr ) ,
#'                  nest = TRUE )
#'
#' # subset
#' flowdes <- subset( flowdes , v2007 == 1 )
#'
#' # means
#' svymean( ~factor( vd4002 ) , flowdes , na.rm = TRUE ) # drops observations with missing in any rounds
#'
#' @export
subset.survflow.design <-
  function (x , subset , rounds = 0 , ... ) {

    if ( is.null( rounds ) ) rounds <- seq_along( x$variables ) - 1
    if ( is.numeric( rounds ) ) {
      rounds <- unique( as.integer( rounds) )
      if ( !all( rounds %in% c( seq_along( x$variables ) - 1 ) ) ) stop( "invalid index; check ?subset.survflow.design for examples.")
    } else stop( "invalid index; check ?update.survflow.design for examples.")

    if ( "svyrep.design" %in% class(x) ) {
      e <- substitute( subset )
      r <- rep( TRUE , nrow( x$variables[[1]] ) )
      for ( z in rounds ) {
        s <- eval(e, x$variables[[z+1]], parent.frame())
        r <- r & s & !is.na(s)
      }
      pwt <- x$pweights
      if (is.data.frame(pwt)) pwt <- pwt[[1]]
      x$pweights <- pwt[r]
      x$repweights <- x$repweights[r, , drop = FALSE]
      if (!is.null(x$selfrep))
        x$selfrep <- x$selfrep[r]
      x$variables <- lapply( x$variables , function(zz) zz[r, , drop = FALSE] )
      x$degf <- NULL
      x$degf <- survey::degf(x)
    } else {
      e <- substitute( subset )
      r <- rep( TRUE , nrow( x$variables[[1]] ) )
      for ( z in rounds ) {
        s <- eval(e, x$variables[[z+1]], parent.frame())
        r <- r & s & !is.na(s)
      }
      x$prob[!r] <- Inf
    }
    return( x )

  }

