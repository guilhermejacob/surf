#' Subset of survflow
#'
#' Restrict a survey flow design to a subpopulation, keeping the original survey flow design information about number of clusters, strata.
#'
#' @method subset survflow.design
#'
#' @param x  a survflow design object
#' @param subset	 An expression specifying the subpopulation
#' @param which.round  a vector of integers indicating which round to apply subsetting condition. Defaults to \code{which.round = 0}, and the expression is applied to the first round.
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
  function (x , subset , which.round = 0 , ... ) {

    if ( is.null( which.round ) ) which.round <- unique( c( 0 , seq_along( attr( x , "data.pairs") ) ) )

    if ( is.numeric( which.round ) ) {

      which.round <- unique( as.integer( which.round) )
      if ( !all( which.round %in% c( 0 , seq_along( attr( x , "data.pairs") ) ) ) ) stop( "invalid index; check ?update.survflow.design for examples.")

      for (z in which.round ) {

        e <- substitute(subset)
        if ( z == 0 ) {
          r <- eval(e, x$variables, parent.frame())
        } else {
          r <- eval(e, attr(x, "data.pairs")[[z]], parent.frame())
        }
        r <- r & !is.na(r)
        x <- x[r, ]
      }
    } else stop( "invalid index; check ?update.survflow.design for examples.")

    return( x )

  }

