#' Flow estimation using complex surveys
#'
#' Specify a complex survey design for flow estimation.
#'
#' @param ids 	Formula or data frame specifying cluster ids from largest level to smallest level, ~0 or ~1 is a formula for no clusters.
#' @param probs	 Formula or data frame specifying cluster sampling probabilities
#' @param strata  Formula or vector specifying strata, use NULL for no strata
#' @param variables	 Formula or data frame specifying the variables measured in the survey. If NULL, the data argument is used.
#' @param fpc  Finite population correction: see Details below
#' @param weights	 Formula or vector specifying sampling weights as an alternative to prob. Notice: longitudinal weights.
#' @param data.list	 list of data frames to look up variables in the formula arguments.
#' @param nest  If TRUE, relabel cluster ids to enforce nesting within strata
#' @param check.strata	 If TRUE, check that clusters are nested in strata
#'
#' @details The arguments of this function are those of \code{\link[survey]{svydesign}},
#' except for the \code{data.list}, which is a list of \code{data.frames}.
#'
#' @return Object of class "\code{survflow.design}", which are \link[survey]{svydesign} based on the initial survey
#' with a \code{data.pairs} attribute containing the data for subsequent surveys.
#'
#' @author Guilherme Jacob
#'
#' @seealso \code{\link[survey]{svydesign}}
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
#' @references Gutierrez, A., Trujillo, L. \& Silva, N. (2014). The estimation of gross ows in complex surveys with random nonresponse,
#' Survey Methodology 40(2), pp. 285-321.
#'
#' Lumley, Thomas S. (2010). Complex Surveys: A Guide to Analysis Using R. Wiley Publishing.
#'
#' @keywords survey
#'
#' @export
svyflowdesign <- function(ids, probs = NULL, strata = NULL, variables = NULL,
                          fpc = NULL, data.list = NULL, nest = FALSE, check.strata = !nest , weights = NULL ){

  # test data input
  if ( length( data.list ) < 2 ) stop( "`data.list` argument must have at least 2 datasets." )
  if ( !all( sapply( data.list , class) %in% "data.frame" ) ) stop( "`data.list` argument must be a list of data.frames" )

  # create design on the first dataset
  res <- survey::svydesign( ids, probs = probs, strata = strata, variables = variables,
                    fpc = fpc, data= data.list[[1]], nest = nest, check.strata = !nest , weights = weights )

  # correct call
  res$call <- match.call()

  # add data pairs
  attr( res , "data.pairs" ) <- data.list[-1]

  # change class
  class( res ) <- c( "survflow.design" , class( res ) )

  # return object
  return( res )

}
