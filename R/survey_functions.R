#' @importFrom stats model.frame model.matrix na.pass printCoefmat terms terms.formula update weights formula as.formula
#' @importFrom methods is

#' @export
as.surfrepdesign <- function( design , ... ) {
  design <- survey::as.svrepdesign( design = design , ... )
  class( design ) <- c( "survflow.design" , class( design ) )
  design
}
