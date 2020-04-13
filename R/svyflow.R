#' @name svyflow
#' @title Gross flow estimation between categories
#'
#' @description Compute gross flows for data from complex surveys with repeated samples.
#'
#' @param x  a formula indicating a \emph{factor} variable.
#' @param design  surflow.design object
#' @param rounds  a vector of integers indicating which round to use. Defaults to \code{rounds = c(0,1)}.
#' @param model  model for non-response. Possibilities: \code{"A", "B", "C", "D"}. Defaults to \code{model = "A"}.
#' @param ...  future expansion.
#'
#' @details The \code{na.rm} option should be used cautiously. Usually, \code{NA} encoding has two possible meanings: a \emph{missing} information or
#' a \emph{not applicable} information. If the former should be dropped, while the latter should not. By default, if the table has an \code{NA}, it will
#' return an table with \code{NAs}. If \code{na.rm = TRUE}, the \code{NA} encoded responses are assumed to be missing at random, and Rojas et al. (2014)
#' method is applied.
#'
#' It is important to distinguish missing responses from unnaplicable responses. This is feasible by  subsetting the design
#' for only applicable responses (with actual missing responses, if that is the case). For instance, suppose that we have two encoded variables:
#' (a) employed/unemployed; and (b) type of contract, with NAs if the response is missing or is unnaplicable. The answers of (a) define what are
#' the applicable responses for responses for (b). In this case, if we are going to analyze contract changes from people employed in the first round,
#' we should filter for people emplyed in both the first \emph{and} second rounds. This can be done using \code{subset}. Then, the remaining \code{NAs} are
#' actual missing responses.
#'
#' @return Objects of class "flowstat", which are tables with a "var" attribute giving the variance and a "statistic" attribute giving the type of flow.
#'
#' These objects have methods for coef, vcov, SE, and cv.
#'
#' @author Guilherme Jacob
#'
#' @examples
#'
#' @references ROJAS, H. A. G.; TRUJILLO, L.; SILVA, P. L. N. The estimation of gross flows in complex surveys with random nonresponse.
#' \emph{Survey Methodology}, v. 40, n. 2, p. 285–321, dec. 2014. URL \url{https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X201400214113}.
#'
#' LUMLEY, T. \emph{Complex Surveys:} A guide to analysis using R.
#' Hoboken: John Wiley & Sons, 2010. (Wiley Series in Survey Methodology). ISBN 978-0-470-28430-8.
#'
#' @examples
#' # load data
#' data( "artificial" )
#'
#' # create surf design object
#' flowdes <-
#'   sfydesign( ids = ~0 ,
#'                  probs = ~ prob ,
#'                  data = list( df0 , df1 ) ,
#'                  nest = TRUE )
#'
#' # flow estimates
#' estflows <- svyflow( ~y , design = flowdes )
#' coef( estflows$muij )
#' SE( estflows$muij )
#'
#' @export
#' @rdname svyflow
#' @method svyflow survey.design2
svyflow.survey.design2 <- function( x , design , rounds , model , ... ){

  # Collect sample data and put in single data frame
  xx <- lapply( design$variables[ rounds + 1 ] , function( z ) stats::model.frame( x , data = z , na.action = stats::na.pass ) )
  xx <- do.call( cbind , xx )

  # Test column format
  if ( !all( sapply( xx , is.factor ) ) ) stop( "this function is only valid for factors." )
  # Gets levels of factors for both time periods
  xlevels <- lapply( xx , function(zz) levels( zz ) )
  # Gets dimension of variable for which flows are to be estimated
  ll <- lapply( xlevels , function(zz) length( zz ) )
  if ( !all( ll[[1]] == ll[[2]] && xlevels[[1]] == xlevels[[2]] ) ) stop( "inconsistent categories across rounds." )
  # When levels are the same across periods, returns unique levels of variable for which flows are to be estimated
  xlevels <- unique( unlist(xlevels) )
  # Gets dimension of variable for which flows are to be estimated
  ll <- unique( unlist(ll) )
  has.order <- ifelse( all( sapply( xx , is.ordered ) ) , TRUE , FALSE )

  # Revise names of columns on sample response data set
  colnames( xx ) <- paste0( "round" , seq_along( colnames( xx ) ) - 1 , ":" , colnames( xx ) )

  # Collect weights
  ww <-  stats::weights( design )

  # model fitting
  mfit <- surf:::ipf( xx , ww , model = model )

  # variance estimation
  mvar <- ipf_variance( xx , ww , res = mfit , design = design )

  # build results list
  res <- sapply( c( "psi" , "rhoRR" , "rhoMM" , "eta" , "pij" , "muij" ) , function(z) {
    if ( z %in% c( "psi" , "rhoRR" , "rhoMM" , "eta" ) ) {
      this_stats <- mfit[[z]]
      attr( this_stats , "var" ) <- mvar[[z]]
      names( attr( this_stats , "var" ) ) <- if ( length( attr( this_stats , "var" ) ) > 1 ) xlevels else z
      class( this_stats ) <- "svystat"
      attr( this_stats , "statistic" ) <- z
    } else if ( z %in% c( "pij" , "muij" ) ) {
      this_stats <- mfit[[z]]
      attr( this_stats , "var" ) <- mfit[[z]]
      attr( this_stats , "var" )[,] <- mvar[[z]][,]
      class( this_stats ) <- "svymstat"
      attr( this_stats , "statistic" ) <- z
    }
    return(this_stats)
  } , simplify = FALSE )

  # create final object
  rval <- res
  rval$model <- mfit$model
  class(rval) <- "flowstat"
  attr( rval , "rounds" )    <- rounds
  attr( rval , "formula" )   <- x
  attr( rval , "has.order" )   <- has.order
  rval

}

#' @export
#' @rdname svyflow
#' @method svyflow svyrep.design
svyflow.svyrep.design <- function( x , design , rounds , model , ... ){


  # Collect sample data and put in single data frame
  xx <- lapply( design$variables[ rounds + 1 ] , function( z ) stats::model.frame( x , data = z , na.action = stats::na.pass ) )
  xx <- do.call( cbind , xx )

  # Test column format
  if ( !all( sapply( xx , is.factor ) ) ) stop( "this function is only valid for factors." )
  # Gets levels of factors for both time periods
  xlevels <- lapply( xx , function(zz) levels( zz ) )
  # Gets dimension of variable for which flows are to be estimated
  ll <- lapply( xlevels , function(zz) length( zz ) )
  if ( !all( ll[[1]] == ll[[2]] && xlevels[[1]] == xlevels[[2]] ) ) stop( "inconsistent categories across rounds." )
  # When levels are the same across periods, returns unique levels of variable for which flows are to be estimated
  xlevels <- unique( unlist(xlevels) )
  # Gets dimension of variable for which flows are to be estimated
  ll <- unique( unlist(ll) )
  has.order <- ifelse( all( sapply( xx , is.ordered ) ) , TRUE , FALSE )

  # Revise names of columns on sample response data set
  colnames( xx ) <- paste0( "round" , seq_along( colnames( xx ) ) - 1 , ":" , colnames( xx ) )

  # Collect weights
  ww <-  stats::weights( design , "sampling" )

  # model fitting
  mfit <- surf:::ipf( xx , ww , model = model )

  # get replication weights
  wr <- stats::weights( design , "analysis" )

  # calculate replicates
  lres <- lapply( seq_len(ncol(wr)) , function( irep , model = mfit$model ) {
    surf:::ipf( xx , wr[,irep] , model = model , starting.values = mfit )
  } )

  # collect replicates
  res1 <- sapply( c( "psi" , "rhoRR" , "rhoMM" , "eta" ) , function(z) {
    qq <- lapply( lres , `[[` , z )
    qq <- do.call( rbind , qq )
    vmat <- survey::svrVar( qq , design$scale , design$rscales , mse = design$mse , coef = matrix( t( mfit[[z]] ) ) )
    diag( vmat )
  } , simplify = FALSE )
  res2 <- sapply( c( "pij" , "muij" ) , function(z) {
    qq <- lapply( lres , `[[` , z )
    qq <- abind::abind( qq , along = 3 )
    nr <- dim(qq)[1]
    qq <- lapply( seq_len(nr) , function(zz) t( qq[zz,,] ) )
    qq <- do.call( cbind , qq )
    vmat <- survey::svrVar( qq , design$scale , design$rscales , mse = design$mse , coef = matrix( t( mfit[[z]] ) ) )
    vmat <- diag( vmat )
    matrix( vmat , byrow = TRUE , nrow = nr )
  } , simplify = FALSE )
  mvar <- c( res1 , res2 )

  # build results list
  res <- sapply( c( "psi" , "rhoRR" , "rhoMM" , "eta" , "pij" , "muij" ) , function(z) {
    if ( z %in% c( "psi" , "rhoRR" , "rhoMM" , "eta" ) ) {
      this_stats <- mfit[[z]]
      attr( this_stats , "var" ) <- mvar[[z]]
      names( attr( this_stats , "var" ) ) <- if ( length( attr( this_stats , "var" ) ) > 1 ) xlevels else z
      class( this_stats ) <- "svystat"
      attr( this_stats , "statistic" ) <- z
    } else if ( z %in% c( "pij" , "muij" ) ) {
      this_stats <- mfit[[z]]
      attr( this_stats , "var" ) <- mfit[[z]]
      attr( this_stats , "var" )[,] <- mvar[[z]][,]
      class( this_stats ) <- "svymstat"
      attr( this_stats , "statistic" ) <- z
    }
    return(this_stats)
  } , simplify = FALSE )

  # create final object
  rval <- res
  rval$model <- mfit$model
  class(rval) <- "flowstat"
  attr( rval , "rounds" )    <- rounds
  attr( rval , "formula" )   <- x
  attr( rval , "has.order" )   <- has.order
  rval

}

#' @export
#' @rdname svyflow
#' @method svyflow surflow.design
svyflow.surflow.design <- function( x , design , rounds = c(0,1) , model = c("A","B","C","D") , ... ) {
  model <- match.arg( model , several.ok = FALSE )
  if ( !all( rounds %in% c(0, seq_along( design$variables ) ) ) ) stop( "rounds not in range." )
  NextMethod( "svyflow" , design = design , rounds = rounds , model = model , ... )
}

#' @export
svyflow <- function( x , design , ... ) {
  # test valid arguments
  if ( class(x) != "formula" ) stop( "x must be a formula." )
  if ( length(x) != 2 ) stop( "x must be a single variable." )
  UseMethod( "svyflow" , design )
}
