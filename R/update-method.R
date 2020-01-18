#' Add variables to a survflow design
#'
#' Update the data variables in a survflow design, either with a formula for a new set of variables or with an expression for variables to be added.
#'
#' @method update survflow.design
#'
#' @param object  a survflow design object
#' @param ...	 Arguments tag=expr add a new variable tag computed by evaluating expr in the survey data
#' @param rounds  a vector of integers indicating which dataset to be updated. Defaults to \code{rounds = NULL}, and the expression is applied to all datasets.
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
#' # update
#' flowdes <- update( flowdes , id_ocup = ifelse( vd4002 ==1 , "ocup" , "desocup" ) )
#'
#' # means
#' svymean( ~factor( id_ocup ) , flowdes , na.rm = TRUE ) # drops observations with missing in any rounds
#'
#' @export
update.survflow.design <-
  function ( object, ... , rounds = NULL ) {

    if ( is.null( rounds ) ) rounds <- unique( c( 0 , seq_along( attr( object , "data.pairs") ) ) )

    if ( is.numeric( rounds ) ) {

      rounds <- unique( as.integer( rounds) )
      if ( !all( rounds %in% c( seq_along( design$variables ) - 1 ) ) ) stop( "invalid index; check ?update.survflow.design for examples.")

      for (z in rounds ) {

        dots <- substitute(list(...))[-1]
        newnames <- names(dots)

        for (j in seq(along = dots)) { object$variables[[z+1]][, newnames[j]] <- eval(dots[[j]], object$variables[[z+1]] , parent.frame()) }

      }
    } else stop( "invalid index; check ?update.survflow.design for examples.")

    return( object )

  }

