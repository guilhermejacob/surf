#' Flow estimation using complex surveys
#'
#' Specify a complex survey design for flow estimation.
#'
#' @param ids 	Formula or data frame specifying cluster ids from largest level to smallest level, ~0 or ~1 is a formula for no clusters.
#' @param probs	 Formula or data frame specifying cluster sampling probabilities
#' @param strata  Formula or vector specifying strata, use NULL for no strata
#' @param fpc  Finite population correction: see Details in \link[survey]{svydesign}.
#' @param weights	 Formula or vector specifying sampling weights as an alternative to prob. Notice: longitudinal weights.
#' @param data.list	 list of data frames to look up variables in the formula arguments.
#' @param nest  If TRUE, relabel cluster ids to enforce nesting within strata
#' @param check.strata	 If TRUE, check that clusters are nested in strata
#'
#' @details The arguments of this function are those of \code{\link[survey]{svydesign}},
#' except for the \code{data.list}, which is a list of \code{data.frames}.
#'
#' The first \code{data.frame} in \code{data.list} must contain the survey design information columns; i.e., clusters, strata, sampling probabilities, etc.
#'
#' @return Object of class \code{surflow.design}, which are \link[survey]{svydesign} with
#' a \code{data.list} attribute containing the data for each survey round.
#'
#' @author Guilherme Jacob
#'
#' @seealso \code{\link[survey]{svydesign}}
#'
#' @examples
#' # load data
#' data( "artificial" )
#'
#' # create surflow design object
#' flowdes <-
#' sfydesign( ids = ~ 1 ,
#'                probs = ~ prob ,
#'                data.list = list( dfa0 , dfa1 ) ,
#'                nest = TRUE )
#'
#' # describe object
#' summary( flowdes )
#'
#' @references Gutierrez, A., Trujillo, L. \& Silva, N. (2014). The estimation of gross flows in complex surveys with random nonresponse,
#' Survey Methodology 40(2), pp. 285-321.
#'
#' Lumley, Thomas S. (2010). Complex Surveys: A Guide to Analysis Using R. Wiley Publishing.
#'
#' @keywords survey
#'
#' @export
sfydesign <- function(ids, probs = NULL, strata = NULL, fpc = NULL, data.list = NULL, nest = FALSE, check.strata = !nest , weights = NULL ){

  # test input
  if ( length( data.list ) < 2 ) stop( "data.list argument must have at least 2 datasets." )
  if ( !all( sapply( data.list , class) %in% "data.frame" ) ) stop( "data.list argument must be a list of data.frames" )
  if ( length( unique( sapply( data.list , nrow ) ) ) > 1 ) stop( "number of observations varies across data.frames. check pairing." )

  # create design on the first dataset
  res <- survey::svydesign( ids, probs = probs, strata = strata, variables = NULL,
                            fpc = fpc, data= data.list[[1]], nest = nest, check.strata = !nest , weights = weights )

  # correct call
  res$call <- match.call()

  # adjusts variables
  res$variables <- data.list

  # change class
  class( res ) <- c( "surflow.design" , class( res ) )

  # return object
  return( res )

}
