#' Add variables to a survflow design
#'
#' Update the data variables in a survflow design, either with a formula for a new set of variables or with an expression for variables to be added.
#'
#' @method update survflow.design
#'
#' @param object  a survflow design object
#' @param which.round  a vector of integers indicating which dataset to be updated. Defaults to \code{which.round = NULL}, and the expression is applied to all datasets.
#' @param ...	 Arguments tag=expr add a new variable tag computed by evaluating expr in the survey data
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
  function (object, which.round = NULL , ... ) {

    if ( is.null( which.round ) ) which.round <- unique( c( 0 , seq_along( attr( object , "data.pairs") ) ) )

    if ( is.numeric( which.round ) ) {

      which.round <- unique( as.integer( which.round) )
      if ( !all( which.round %in% c( 0 , seq_along( attr( object , "data.pairs") ) ) ) ) stop( "invalid index; check ?update.survflow.design for examples.")

      for (z in which.round ) {

        dots <- substitute(list(...))[-1]
        newnames <- names(dots)

        if ( z == 0 ) {
          for (j in seq(along = dots)) { object$variables[, newnames[j]] <- eval(dots[[j]], object$variables, parent.frame()) }
        } else {
          for (j in seq(along = dots)) { attr( object, "data.pairs" )[[z]][, newnames[j]] <- eval(dots[[j]], attr( object, "data.pairs" )[[z]], parent.frame()) }
        }
      }
    } else stop( "invalid index; check ?update.survflow.design for examples.")

    return( object )

  }

