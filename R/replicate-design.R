#' Convert surflow design to use replicate weights
#'
#' Creates a replicate-weights survey flow design object from a traditional strata/cluster survey flow design object.
#' JK1 and JKn are jackknife methods, BRR is Balanced Repeated Replicates and Fay is Fay's modification of this,
#' bootstrap is Canty and Davison's bootstrap, subbootstrap is Rao and Wu's (n-1) bootstrap,
#' and mrbbootstrap is Preston's multistage rescaled bootstrap.
#'
#' @param design 	Object of class surflow.design
#' @param ...	 Arguments passed to \code{\link[survey]{as.svrepdesign}}
#'
#' @details The arguments of this function are those of \code{\link[survey]{as.svrepdesign}},
#'
#' @return Object of class \code{surflow.design}, \link[survey]{svrepdesign} with
#' a \code{data} attribute containing the data for each survey round.
#'
#' @author Guilherme Jacob
#'
#' @seealso \code{\link[survey]{as.svrepdesign}}
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
#' # transform in replicate design
#' flowdes_rep <- as.surfrdesign( flowdes , type = "bootstrap" , replicates = 50 )
#'
#' @references LUMLEY, T. \emph{Complex Surveys:} A guide to analysis using R.
#' Hoboken: John Wiley & Sons, 2010. (Wiley Series in Survey Methodology). ISBN 978-0-470-28430-8.
#'
#' @keywords survey
#'
#' @export
#' @export
as.surfrdesign <- function( design , ... ) {
  design <- survey::as.svrepdesign( design = design , ... )
  class( design ) <- c( "surflow.design" , class( design ) )
  attr( design , "fullpweights" ) <- design$pweights
  attr(design , "fullrepweights" ) <- design$repweights
  design
}
