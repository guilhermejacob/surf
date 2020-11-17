#' @name svyflow
#' @title Gross flow estimation between categories
#'
#' @description Compute gross flows for data from complex surveys with repeated samples.
#'
#' @param x  a formula indicating a \emph{factor} variable.
#' @param design  surflow.design object
#' @param model  model for non-response. Possibilities: \code{"A", "B", "C", "D"}. Defaults to \code{model = "A"}.
#' @param tol  Tolerance for iterative proportional fitting. Defaults to \code{1e-4}.
#' @param maxit  Maximum number of iterations for iterative proportional fitting. Defaults to \code{1000}.
#' @param verbose  Print proportional fitting iterations. Defaults to \code{verbose = FALSE}.
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
#' we should filter for people employed in both the first \emph{and} second rounds. This can be done using \code{subset}. Then, the remaining \code{NA} are
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
#'
#' # load library
#' library( survey )
#' library( surf )
#'
#' # load data
#' data( "LFS79.0809" )
#'
#' # create surf design object
#' lfs.des <- svydesign( ids = ~0 , probs = ~ prob , data = LFS79.0809 , nest = TRUE )
#'
#' # flow estimates
#' estflows <- svyflow( ~y1+y2 , design = lfs.des )
#' coef( estflows$muij )
#' SE( estflows$muij )
#'
#' @export
#' @rdname svyflow
#' @method svyflow survey.design2
svyflow.survey.design2 <- function( x , design , model = c("A","B","C","D") , tol = 1e-4 , maxit = 1000 , verbose = FALSE , ... ){

  # test values
  model <- match.arg( model , several.ok = FALSE )

  # collect sample data and put in single data frame
  xx <- stats::model.frame( x , data = design$variables , na.action = stats::na.pass )

  # test column format
  if ( !all( sapply( xx , is.factor ) ) ) stop( "this function is only valid for factors." )

  # Gets levels of factors for both time periods
  xlevels <- lapply( xx , function(zz) levels( zz ) )

  # gets dimension of variable for which flows are to be estimated
  ll <- lapply( xlevels , function(zz) length( zz ) )
  if ( !all( ll[[1]] == ll[[2]] & xlevels[[1]] == xlevels[[2]] ) ) stop( "inconsistent categories across rounds." )

  # when levels are the same across periods, returns unique levels of variable for which flows are to be estimated
  xlevels <- unique( unlist(xlevels) )

  # check for ordered categories
  ll <- unique( unlist(ll) )
  has.order <- ifelse( all( sapply( xx , is.ordered ) ) , TRUE , FALSE )

  # collect weights
  ww <-  stats::weights( design )

  # estimate counts
  Amat <- stats::xtabs( c(ww,0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE , drop.unused.levels = FALSE )

  # non-response cells
  Nij <- Amat[ -nrow( Amat ) , -ncol( Amat ) ]
  Ri <- Amat[ , ncol(Amat) ][ - nrow( Amat ) ]
  Cj <- Amat[ nrow( Amat ) , ][ - ncol( Amat ) ]
  M <- Amat[ nrow( Amat ) , ncol( Amat ) ]

  # treat full response
  if ( all( c( Ri , Cj , M ) == 0 ) & any( Nij > 0 ) ) {

    # issue warning
    warning( "counts show full response. ignoring model.")

    # remove borders from table
    Amat <- Amat[ -nrow( Amat ) , -ncol( Amat ) ]

    # model fitting
    mfit <- frf( Amat )
    mfit$delta <- mfit$gamma - mfit$eta

    # variance estimation
    mvar <- frf_variance( xx , ww , res = mfit , design = design )

    # build results list
    res <- sapply( c( "eta" , "pij" , "muij" , "gamma" , "delta" ) , function(z) {
      if ( z %in% c( "psi" , "rho" , "tau" , "eta" , "gamma" , "delta" ) ) {
        this_stats <- mfit[[z]]
        attr( this_stats , "var" ) <- mvar[[z]]
        names( this_stats ) <- if ( length( this_stats ) > 1 ) xlevels else z
        colnames( attr( this_stats , "var" ) ) <- rownames( attr( this_stats , "var" ) ) <- if ( length( attr( this_stats , "var" ) ) > 1 ) xlevels else z
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
    rval <- res[ c( "psi" , "rho" , "tau" , "eta" , "gamma" , "pij" , "muij" , "delta" ) ]
    rval$model <- "Full Response"
    class(rval) <- "flowstat"
    attr( rval , "formula" )   <- x
    attr( rval , "has.order" )   <- has.order
    attr( rval , "iter" )   <- NA
    attr( rval , "observed.counts" )   <- mfit$observed.counts

    # return final object
    return( rval )

  }

  # test for zero counts
  if ( any( Amat <= 0 ) ) {
    # issue warning
    warning( "stopping. some cells had zero counts. consider collapsing categories.")
    return( Amat )
  }

  # model fitting
  mfit <- ipf( Amat , model = model , tol = tol , maxit = maxit , verbose = verbose )

  # variance estimation
  mvar <- ipf_variance( xx , ww , res = mfit , design = design , rp.variance = TRUE )

  # build results list
  res <- sapply( c( "psi" , "rho" , "tau" , "eta" , "pij" , "muij" , "gamma" , "delta" ) , function(z) {
    if ( z %in% c( "psi" , "rho" , "tau" , "eta" , "gamma" , "delta" ) ) {
      this_stats <- mfit[[z]]
      attr( this_stats , "var" ) <- mvar[[z]]
      names( this_stats ) <- if ( length( this_stats ) > 1 ) xlevels else z
      colnames( attr( this_stats , "var" ) ) <- rownames( attr( this_stats , "var" ) ) <- if ( length( attr( this_stats , "var" ) ) > 1 ) xlevels else z
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
  rval <- res[ c( "psi" , "rho" , "tau" , "eta" , "gamma" , "pij" , "muij" , "delta" ) ]
  rval$model <- mfit$model
  class(rval) <- "flowstat"
  attr( rval , "formula" )   <- x
  attr( rval , "has.order" )   <- has.order
  attr( rval , "iter" )   <- mfit$iter
  attr( rval , "unadj.chisq" )   <- mfit$unadj.chisq
  attr( rval , "adj.chisq" )   <- mvar$adj.chisq
  attr( rval , "observed.counts" )   <- mfit$observed.counts

  # return final object
  rval

}

#' @export
svyflow <- function( x , design , model = c( "A","B","C","D") , ... ) {

  # test values
  model <- match.arg( model , several.ok = FALSE )

  # test valid arguments
  if ( class( x ) != "formula" ) stop( "x must be a formula." )
  if ( length( as.character( x ) ) != 2 ) stop( "x must be a one-sided formula." )
  if ( length( ( strsplit( as.character( x )[[2]] , " \\+ " ) )[[1]] ) != 2 ) stop("only two-way tables at the moment.")
  if ( ncol( attr( terms( x ) , "factors" ) ) != 2 ) stop("only two-way tables at the moment.")
  UseMethod( "svyflow" , design )

}
