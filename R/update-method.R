#' @title Add variables to a surflow design
#'
#' @description Update the data variables in a surflow design, either with a formula for a new set of variables or with an expression for variables to be added.
#'
#' @param object  a surflow design object
#' @param ...	 Arguments tag=expr add a new variable tag computed by evaluating expr in the survey data
#' @param rounds  a vector of integers indicating which dataset to be updated. Defaults to \code{rounds = NULL}, and the expression is applied to all datasets.
#'
#' @examples
#' # load data
#' data( "artificial" )
#'
#' # create surflow design object
#' flowdes <-
#' sfydesign( ids = ~ 1 ,
#'                probs = ~ prob ,
#'                data = list( dfa0 , dfa1 ) ,
#'                nest = TRUE )
#'
#' # update
#' flowdes <-
#'      update( flowdes ,
#'              idstatus = factor( v0 ,
#'                                 levels = 1:3 ,
#'                                 labels = c( "A" , "B" , "C" ) ) )
#'
#' # evaluate function
#' svyflow( ~idstatus , flowdes )
#'
#' @method update surflow.design
#' @export
update.surflow.design <-
  function ( object, ... , rounds = NULL ) {

    if ( is.null( rounds ) ) rounds <- seq_along( object$variables ) - 1

    if ( is.numeric( rounds ) ) {

      rounds <- unique( as.integer( rounds) )
      if ( !all( rounds %in% c( seq_along( object$variables ) - 1 ) ) ) stop( "invalid index; check ?update.surflow.design for examples.")

      for (z in rounds ) {

        dots <- substitute(list(...))[-1]
        newnames <- names(dots)
        for (j in seq(along = dots)) { object$variables[[z+1]][, newnames[j]] <- eval(dots[[j]], object$variables[[z+1]] , parent.frame()) }

      }
    } else stop( "invalid index; check ?update.surflow.design for examples.")

    return( object )

  }

