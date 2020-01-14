#' @name surveysummary
#' @title Summary statistics for repeated sample surveys
#'
#' @description Compute means and totals for data from complex surveys with repeated samples.
#'
#' @aliases svytotal svymean
#'
#' @param x  A formula
#' @param design  survflow.design object
#' @param which.round  a vector of integers indicating which round to use. Defaults to \code{which.round = NULL} and calculation is repeated in all rounds.
#' @param ...  arguments to be passed to \code{\link[survey]{surveysummary}} function.
#'
#' @details The arguments of this function are those of \code{\link[survey]{surveysummary}}.
#'
#' @return Objects of class "svystat" or "svrepstat", which are vectors with a "var" attribute giving the variance and a "statistic" attribute giving the name of the statistic.
#'
#' These objects have methods for vcov, SE, coef, confint, svycontrast.
#'
#' @author Guilherme Jacob
#'
#' @seealso \code{\link[survey]{surveysummary}}
#'
#' @examples
#' # load data
#' data( "initial" )
#' data( "final" )
#'
#' # create survflow design object
#' flowdes <-
#' svyflowdesign( ids = ~ upa ,
#'                strata = ~ estrato ,
#'                probs = ~ longprob ,
#'                data.list = list( f.qtr , s.qtr ) ,
#'                nest = TRUE )
#'
#' # describe object
#' summary( flowdes )
#'
#' @references Lumley, Thomas S. (2010). Complex Surveys: A Guide to Analysis Using R. Wiley Publishing.
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
#' # totals
#' svytotal( ~factor( vd4002 ) , flowdes , na.rm = TRUE ) # drops observations with missing in any rounds
#' svytotal( ~factor( vd4002 ) , flowdes , na.rm = TRUE , which.dataset = 0 ) # drops observations with missing in the first round
#' svytotal( ~factor( vd4002 ) , flowdes , na.rm = TRUE , which.dataset = 1 ) # drops observations with missing in the second round
#'
#' # means
#' svymean( ~factor( vd4002 ) , flowdes , na.rm = TRUE ) # drops observations with missing in any rounds
#' svymean( ~factor( vd4002 ) , flowdes , na.rm = TRUE , which.dataset = 0 ) # drops observations with missing in the first round
#' svymean( ~factor( vd4002 ) , flowdes , na.rm = TRUE , which.dataset = 1 ) # drops observations with missing in the second round
#'
#' @export
#' @rdname surveysummary
#' @method svytotal survflow.design
svytotal.survflow.design <- function( x , design , which.round = NULL , ... ){

  # collect data
  x <- lapply( c( list( design$variables ) , attr( design , "data.pairs" ) ) , function( z ) stats::model.frame( x , data = z , na.action = na.pass ) )
  x <- do.call( cbind , x )

  # add time frame
  colnames( x ) <- paste0( "round" , seq_along( colnames( x ) ) , ":" , colnames( x ) )

  # filter rounds
  if ( is.null( which.round ) ) { which.round <- unique( c( 0 , seq_along( attr( design , "data.pairs") ) ) ) + 1 } else {
    if ( !is.numeric( which.round ) ) stop( "invalid index; check ?update.survflow.design for examples.")
    which.round <- unique( as.integer( which.round) )
    if ( !all( which.round %in% c( 0 , seq_along( attr( design , "data.pairs") ) ) ) ) stop( "invalid index; check ?update.survflow.design for examples.")
    which.round <- which.round + 1
  }

  # apply function
  class( design ) <- class( design )[-1]
  UseMethod( "svytotal" , object = design , x , `...` = `...` )

}


#' @export
#' @rdname surveysummary
#' @method svymean survflow.design
svymean.survflow.design <- function( x , design , which.round = NULL , ... ){

  # collect data
  x <- lapply( c( list( design$variables ) , attr( design , "data.pairs" ) ) , function( z ) stats::model.frame( x , data = z , na.action = na.pass ) )
  x <- do.call( cbind , x )

  # add time frame
  colnames( x ) <- paste0( "round" , seq_along( colnames( x ) ) , ":" , colnames( x ) )

  # filter rounds
  if ( is.null( which.round ) ) { which.round <- unique( c( 0 , seq_along( attr( design , "data.pairs") ) ) ) + 1 } else {
    if ( !is.numeric( which.round ) ) stop( "invalid index; check ?update.survflow.design for examples.")
    which.round <- unique( as.integer( which.round) )
    if ( !all( which.round %in% c( 0 , seq_along( attr( design , "data.pairs") ) ) ) ) stop( "invalid index; check ?update.survflow.design for examples.")
    which.round <- which.round + 1
  }

  # apply function
  class( design ) <- class( design )[-1]
  UseMethod( "svymean" , object = design , x , `...` = `...` )

}

#' @export
svytotal <- survey::svytotal
#' @export
svymean <- survey::svymean
