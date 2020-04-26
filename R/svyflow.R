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

  # if no missing, use simple contingency table
  if ( !any( is.na( xx[ ww > 0 ,] ) ) ) {

    bigNij <- stats::xtabs( c(ww,0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE , drop.unused.levels = FALSE )
    bigNij <- as.matrix( bigNij[ -nrow( bigNij ) , -ncol( bigNij ) ] )
    eta <- rowSums( bigNij ) / sum( bigNij )
    pij <- bigNij / rowSums( bigNij )
    muij <- sum( ww ) * sweep( pij , 1 , eta , "*" )

    ### variance

    # yy array - see Rojas et al. (2014, p.294)
    yy <- array( 0  , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( xx ) ) )
    for ( r in seq_len( ncol( xx ) ) ) {
      kk <- stats::model.matrix( ~-1+. , data = xx[,r,drop = FALSE] , contrasts.arg = lapply( xx[,r, drop = FALSE ] , stats::contrasts, contrasts=FALSE ) , na.action = stats::na.pass )
      yy[ which( !is.na( xx[ , r ] ) ) , , r ] <- kk ; rm( kk )
    }

    # Create matrix of z variables - see Rojas et al. (2014, p.294)
    zz <- apply( xx , 2 , function(z) as.numeric( !is.na(z) ) )

    # Special variables - see Rojas et al. (2014, p.295)
    vv <- rowSums( yy[,,1] ) * rowSums( yy[,,2] ) + rowSums( yy[,,2] * (1 - zz[,1]) ) + rowSums( yy[,,1] * (1 - zz[,2]) ) + ( 1- zz[,1] ) * ( 1- zz[,2] )
    y1y2 <- array( 0  , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) for ( j in seq_len( ncol( bigNij ) ) ) y1y2[,i,j] <- yy[,i,1] * yy[,j,2]

    # aux stats
    nipij <- sweep( pij , 1 , eta , "*" )

    ### eta

    # Calculate scores for estimating the variance of eta parameters
    u_eta <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_eta[,i] <-
        ( sum( bigNij ) * rowSums( y1y2[,i,] ) - rowSums( bigNij )[i] * rowSums( y1y2[,i,] ) * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) ) / sum( bigNij )^2
    }

    # calculate variance of eta
    eta_var <- survey::svyrecvar( sweep( u_eta , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    eta_var <- diag( eta_var )

    ### pij

    # Calculate scores for estimating the variance of pij parameters
    u_pij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_pij[,i,j] <-
        ( sum( bigNij[i,] ) * y1y2[,i,j] - bigNij[i,j] * rowSums( y1y2[,i, ] ) ) / sum( bigNij[i,] )^2
    }

    # calculate variance of u_pij sum
    lin_pij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_pij[,z,] ) )
    pij_var <- diag( survey::svyrecvar( ww * lin_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    pij_var <- matrix( pij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + sum(ww) * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    # calculate variance of muij
    lin_muij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_muij[,z,] ) )
    muij_var <- diag( survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    muij_var <- matrix( muij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### combine results
    mfit <- list( "psi" = NA , "rhoRR" = NA , "rhoMM" = NA , "eta" = eta , "pij" = pij , "muij" = bigNij )
    mfit$model <- "(No Missing)"
    mvar <- list( "psi" = NA , "rhoRR" = NA , "rhoMM" = NA , "eta" = eta_var , "pij" = pij_var , "muij" = muij_var )

  } else {

    # model fitting
    mfit <- ipf( xx , ww , model = model , tol = list(`...`)[["tol"]] , verbose = list(`...`)[["verbose"]] , starting.values = list(`...`)[["starting.values"]] )

    # variance estimation
    mvar <- ipf_variance( xx , ww , res = mfit , design = design )

  }

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
  attr( rval , "iter" )   <- mfit$iter
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
  mfit <- ipf( xx , ww , model = model , tol = tol , verbose = verbose )

  # get replication weights
  wr <- stats::weights( design , "analysis" )

  # calculate replicates
  lres <- lapply( seq_len(ncol(wr)) , function( irep , model = mfit$model ) {
    ipf( xx , wr[,irep] , model = model , starting.values = mfit )
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
  attr( rval , "iter" )   <- mfit$iter
  rval

}

#' @export
#' @rdname svyflow
#' @method svyflow surflow.design
svyflow.surflow.design <- function( x , design , rounds = c(0,1) , model = c("A","B","C","D") , tol = 1e-6 , verbose = FALSE , starting.values = list( "psi" = NULL , "rhoRR" = NULL , "rhoMM" = NULL, "eta" = NULL , "pij" = NULL ) ) {
  model <- match.arg( model , several.ok = FALSE )
  if ( !all( rounds %in% c(0, seq_along( design$variables ) ) ) ) stop( "rounds not in range." )
  NextMethod( "svyflow" , design = design , rounds = rounds , model = model , tol = tol , verbose = verbose , starting.values = starting.values )
}

#' @export
svyflow <- function( x , design , ... ) {
  # test valid arguments
  if ( class(x) != "formula" ) stop( "x must be a formula." )
  if ( length(x) != 2 ) stop( "x must be a single variable." )
  UseMethod( "svyflow" , design )
}
