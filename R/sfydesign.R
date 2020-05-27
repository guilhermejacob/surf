#' @title Flow estimation using complex surveys
#'
#' @description Specify a complex survey design for flow estimation.
#'
#' @usage sfydesign( ids, probs=NULL, strata=NULL, fpc=NULL, data=NULL, nest=FALSE,
#' check.strata = !nest, weights=NULL, pps=FALSE,...)
#'
#' @param ids 	Formula or data frame specifying cluster ids from largest level to smallest level, ~0 or ~1 is a formula for no clusters.
#' @param probs	 Formula or data frame specifying cluster sampling probabilities
#' @param strata  Formula or vector specifying strata, use NULL for no strata
#' @param fpc  Finite population correction: see Details in \link[survey]{svydesign}.
#' @param weights	 Formula or vector specifying sampling weights as an alternative to prob. Notice: longitudinal weights.
#' @param data	 list of data frames to look up variables in the formula arguments.
#' @param nest  If TRUE, relabel cluster ids to enforce nesting within strata
#' @param check.strata	 If TRUE, check that clusters are nested in strata
#' @param pps	 \code{"brewer"} to use Brewer's approximation for PPS sampling without replacement.
#' \code{"overton"} to use Overton's approximation. An object of class \link[survey]{HR} to use the Hartley-Rao approximation.
#' An object of class \link[survey]{ppsmat} to use the Horvitz-Thompson estimator.
#' @param ...  for future expansion. See \link[survey]{svydesign} for more information.
#'
#' @details The arguments of this function are those of \code{\link[survey]{svydesign}},
#' except for the \code{data}, which, in this case, is a list of \code{data.frames} with \emph{paired observations}; i.e., each
#' row must refer to the same individual across datasets.
#'
#' The first \code{data.frame} in \code{data} must contain the survey design information columns; i.e., clusters, strata, sampling probabilities, etc.
#'
#' @return Object of class \code{surflow.design}, which are \link[survey]{svydesign} with
#' a \code{data} attribute containing the data for each survey round.
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
#'                data = list( df1 , df2 ) ,
#'                nest = TRUE )
#'
#' # describe object
#' summary( flowdes )
#'
#' @references ROJAS, H. A. G.; TRUJILLO, L.; SILVA, P. L. N. The estimation of gross flows in complex surveys with random nonresponse.
#' \emph{Survey Methodology}, v. 40, n. 2, p. 285–321, dec. 2014. URL \url{https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X201400214113}.
#'
#' LUMLEY, T. \emph{Complex Surveys:} A guide to analysis using R.
#' Hoboken: John Wiley & Sons, 2010. (Wiley Series in Survey Methodology). ISBN 978-0-470-28430-8.
#'
#' @keywords survey
#'
#' @export
sfydesign <- function(ids, probs = NULL, strata = NULL, fpc = NULL, data = NULL, nest = FALSE, check.strata = !nest , weights = NULL, pps=FALSE, ... ){

  # test input
  if ( class( data ) != "list" ) stop( "data argument must be a list of data.frames" )
  if ( length( data ) < 2 ) stop( "data argument must have at least 2 datasets." )
  if ( !all( sapply( data , class) %in% "data.frame" ) ) stop( "data argument must be a list of data.frames" )
  if ( length( unique( sapply( data , nrow ) ) ) > 1 ) stop( "number of observations varies across data.frames. check pairing." )

  # preprocess
  data <- lapply( data , function( zz ) {rownames(zz) <- seq_len( nrow(zz)) ; zz} )

  # create design on the first dataset
  res <- survey::svydesign( ids, probs = probs, strata = strata, variables = NULL,
                            fpc = fpc, data= data[[1]], nest = nest, check.strata = !nest , weights = weights , pps = pps , ... )

  # full sampling weights (used to calculate endogenous categories)
  attr( res , "fullprob" ) <- stats::weights( res )

  # correct call
  res$call <- match.call()

  # adjusts variables
  res$variables <- data

  # change class
  class( res ) <- c( "surflow.design" , class( res ) )

  # return object
  return( res )

}
