#' @name svyflow
#' @title Gross flow estimation between (exogenous) categories
#'
#' @description Compute gross flows for data from complex surveys with repeated samples.
#'
#' @param x  A formula indcating a \emph{factor} variable.
#' @param design  surflow.design object
#' @param flow.type  type of flow to estimate: "gross" for counts and "net" for probabilities. Defaults to \code{flow.type = "gross"}.
#' @param rounds  a vector of integers indicating which round to use. Defaults to \code{rounds = c(0,1)}.
#' @param max.iter  number of iterations. Defaults to \code{max.iter = 10}.
#' @param extra  Should initial and transition probabilities be stored? Defaults to \code{extra = FALSE}.
#' @param na.rm  Should missing variables be dropped? Defaults to \code{na.rm = FALSE}. See details for further information.
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
#'                  data = list( dfa0 , dfa1 ) ,
#'                  nest = TRUE )
#'
#' # gross flows
#' gross.flows <- svyflow( ~v0 , design = flowdes , flow.type = "gross" )
#' coef( gross.flows )
#' SE( gross.flows )
#'
#' # net flows
#' net.flows <- svyflow( ~v0 , design = flowdes , flow.type = "net" )
#' coef( net.flows )
#' SE( net.flows )
#'
#' @export
#' @rdname svyflow
#' @method svyflow survey.design2
svyflow.survey.design2 <- function( x , design , flow.type , rounds , max.iter , extra , na.rm , ... ){

  # collect data
  xx <- lapply( design$variables[ rounds + 1 ] , function( z ) stats::model.frame( x , data = z , na.action = stats::na.pass ) )
  xx <- do.call( cbind , xx )

  # test column format
  if ( !all( sapply( xx , is.factor ) ) ) stop( "this function is only valid for factors." )
  xlevels <- lapply( xx , function(zz) levels( zz ) )
  ll <- lapply( xlevels , function(zz) length( zz ) )
  if ( !all( ll[[1]] == ll[[2]] && xlevels[[1]] == xlevels[[2]] ) ) stop( "inconsistent categories across rounds." )
  xlevels <- unique( unlist(xlevels) )
  ll <- unique( unlist(ll) )
  has.order <- ifelse( all( sapply( xx , is.ordered ) ) , TRUE , FALSE )

  # add time frame
  colnames( xx ) <- paste0( "round" , seq_along( colnames( xx ) ) - 1 , ":" , colnames( xx ) )

  # collect weights
  ww <- 1 / design$prob

  # aggregate
  NN <- stats::xtabs( c(ww,0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE , drop.unused.levels = FALSE )
  RR <- NN[ , ncol(NN) ][ - nrow( NN ) ]
  CC <- NN[ nrow(NN) , ][ - ncol( NN ) ]
  MM <- NN[ nrow( NN ) , ncol( NN ) ]
  NN <- as.matrix( NN[ -nrow( NN ) , -ncol( NN ) ] )
  N  <- sum( NN ) + sum( RR ) + sum( CC ) + MM

  # handling missing
  if ( sum( MM , RR , CC ) > 0 && !na.rm ) {
    # format results
    NN[,] <- NA
    rval <- NN
    rval <- as.table(rval)
    class(rval) <- "flowstat"
    attr( rval , "var" )       <- unclass( NN )
    attr( rval , "statistic" ) <- flow.type
    attr( rval , "rounds" )    <- rounds
    attr( rval , "formula" )   <- x
    attr( rval , "has.order" )   <- has.order
    return(rval)
  }

  # model fitting
  res <- ipf( xx , ww )

  # collect results
  eta_i <- res[["eta_i"]]
  p_ij <- res[["p_ij"]]
  N <- res[["N"]]
  nipij <- res[["nipij"]]
  mu_ij <- N * nipij
  if ( any( mu_ij == 0 ) ) warning( "Some flows are equal to zero. Variances are inconsistent." )

  ### variance calculations

  # yy array
  yy <- array( 0  , dim = c( nrow( xx ) , ll , ncol( xx ) ) )
  for ( r in seq_len( ncol( xx ) ) ) {
    kk <- stats::model.matrix( ~-1+. , data = xx[,r,drop = FALSE] , contrasts.arg = lapply( xx[,r, drop = FALSE ] , stats::contrasts, contrasts=FALSE ) , na.action = stats::na.pass )
    yy[ which( !is.na( xx[ , r ] ) ) , , r ] <- kk ; rm( kk )
  }

  # create matrix of z variables
  zz <- apply( xx , 2 , function(z) as.numeric( !is.na(z) ) )

  ### special variables
  vv_k <- rowSums( yy[,,1] ) * rowSums( yy[,,2] ) + rowSums( yy[,,2] * (1 - zz[,1]) ) + rowSums( yy[,,1] * (1 - zz[,2]) ) + ( 1- zz[,1] ) * ( 1- zz[,2] )
  nnij_k <- array( 0  , dim = c( nrow( xx ) , nrow( NN ) , ncol( NN ) ) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) nnij_k[,i,j] <- yy[,i,1] * yy[,j,2]

  ### variance of eta
  ueta_ik <- array( 0 , dim = c( nrow(xx) , nrow(NN) ) )
  for ( i in seq_len( nrow(NN) ) ) {
    ueta_ik[,i] <-
      ( rowSums( nnij_k[,i,] ) + yy[,i,1] * ( 1 - zz[,2] ) ) / eta_i [i] +
      rowSums( sweep( yy[,,2] * ( 1 - zz[,1] ) , 2 , p_ij[i,] / colSums( nipij ) , "*" ) ) +
      (1 - zz[,1]) * (1 - zz[,2])
  }

  jeta_i <- vector( "numeric" , length = nrow(NN) )
  for ( i in seq_len( nrow(NN) ) ) {
    jeta_i[i] <- ( sum( ww * yy[,i,1] * zz[,2] ) - 2 * sum( ww * yy[,i,1] ) ) / ( eta_i[i]^2 )
    - sum( ww * ( 1 - zz[,1] ) * rowSums( sweep( yy[,,2] , 2 , (p_ij[i,] / colSums( nipij ))^2 , "*" ) ) )
  }

  if (extra) {
    # adjusting factor
    lin_eta <- sweep( ueta_ik , 2 , jeta_i , "/" )
    # eta_var <- diag( survey::svyrecvar( ww * lin_eta , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    eta_var <- survey::svyrecvar( ww * lin_eta , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    rm( lin_eta )
    eta_res <- eta_i
    colnames( eta_var ) <- rownames( eta_var ) <-  names( eta_res ) <- xlevels
    class( eta_res ) <- "svystat"
    attr( eta_res , "var" ) <- eta_var
    attr( eta_res , "statistic" ) <- "initial probabilities"
  }

  # variance of pij
  up_ijk <- array( 0 , dim = c( nrow(xx) , nrow(NN) , ncol(NN) ) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    up_ijk[,i,j] <- ( nnij_k[,i,j] / p_ij[i,j] ) +
      ( yy[,i,1] * ( 1 - zz[,2] ) ) +
      ( yy[,j,2] * ( 1 - zz[,1] ) ) * ( eta_i[i] / colSums( nipij )[j] ) +
      ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * eta_i[ i ]
  }

  jp_ij <- array( 0 , dim = c( nrow(NN) , ncol(NN) ) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    jp_ij[i,j] <- - sum( ww * nnij_k[,i,j] ) / ( p_ij[i,j]^2 ) -
      ( ( eta_i[i] / colSums( nipij )[j] )^2 ) *
      sum( ww * yy[,j,2] * ( 1 - zz[,1] ) )
  }

  if (extra) {
    # adjusting factor
    lin_pij <- array( 0 , dim = c( nrow(xx) , nrow(NN) , ncol(NN) ) )
    lin_pij <- sweep( up_ijk , c(2,3) , jp_ij , "/" )
    # for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) lin_pij[ , i,j] <- up_ijk[,i,j] / jp_ij[i,j]
    lin_pij <- do.call( cbind , lapply( seq_len(ncol( NN)) , function( z ) lin_pij[,z,] ) )
    p_var <- diag( survey::svyrecvar( ww * lin_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    p_var <- matrix( p_var , nrow = nrow(NN) , ncol = ncol(NN) , byrow = TRUE )
    rm( lin_pij )
    pij_res <- p_ij
    # colnames( p_var ) <- rownames( p_var ) <-  names( pij_res ) <- xlevels
    class( pij_res ) <- "flowstat"
    attr( pij_res , "var" ) <- p_var
    attr( pij_res , "statistic" ) <- "transition probabilities"
  }
  # # variance of NN
  # lin_nnij <- do.call( cbind , lapply( seq_len(ncol( NN)) , function( z ) nnij_k[,z,] ) )
  # NN_var <- diag( survey::svyrecvar( ww * lin_nnij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
  # NN_var <- matrix( NN_var , nrow = nrow(NN) , ncol = ncol(NN) , byrow = TRUE )
  #
  # # a variables
  # mu_var <- matrix( as.numeric(NA) , nrow = nrow(NN) , ncol = ncol(NN) )
  # for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
  #   a7 <- nipij[i,j]
  #   a8 <- NN[i,j] * p_ij[i,j]
  #   a9 <- NN[i,j] * eta_i[i]
  #   mu_var[i,j] <- (a7^2) * NN_var[i,j] + (a8^2) * eta_var[i] + (a9^2) * p_var[i,j]
  # }

  # alternative strategy
  nipij_k <- array( 0 , dim = c( nrow( xx ) , nrow( NN ) , ncol( NN ) ) )
  muij_k <- array( 0 , dim = c( nrow( xx ) , nrow( NN ) , ncol( NN ) ) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    nipij_k[,i,j] <- p_ij[i,j] * ( ueta_ik[,i] / jeta_i[i] ) + eta_i[i] * ( up_ijk[,i,j] / jp_ij[i,j] )
    muij_k[,i,j] <- nipij[i,j] * vv_k + N * nipij_k[,i,j]
  }
  rm( yy , ueta_ik , jeta_i , up_ijk , jp_ij )
  nipij_k <- do.call( cbind , lapply( seq_len( ncol( NN) ) , function( z ) nipij_k[,z,] ) )
  muij_k <- do.call( cbind , lapply( seq_len( ncol( NN) ) , function( z ) muij_k[,z,] ) )
  nipij_k[ is.na( nipij_k ) | is.infinite( nipij_k ) ] <- 0
  muij_k[ is.na( muij_k ) | is.infinite( nipij_k ) ] <- 0
  if ( flow.type == "gross" ) {
    mu_var <- diag( survey::svyrecvar( ww * muij_k , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
  } else {
    mu_var <- diag( survey::svyrecvar( ww * nipij_k , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
  }
  mu_var <- matrix( mu_var , nrow = nrow( NN ) , ncol = ncol( NN ) , byrow = TRUE )

  # format results
  rval <- if ( flow.type == "gross" ) mu_ij else nipij
  rval <- unclass(rval)
  class(rval) <- "flowstat"
  attr( rval , "var" )       <- unclass(mu_var)
  dimnames( attr( rval , "var" ) ) <- dimnames( mu_ij )
  attr( rval , "statistic" ) <- flow.type
  attr( rval , "rounds" )    <- rounds
  attr( rval , "formula" )   <- x
  attr( rval , "has.order" )   <- has.order
  attr( rval , "eta" ) <- if (extra) eta_res else NULL
  attr( rval , "pij" ) <- if (extra) pij_res else NULL
  rval

}

#' @export
#' @rdname svyflow
#' @method svyflow svyrep.design
svyflow.svyrep.design <- function( x , design , flow.type , rounds , max.iter , extra , na.rm , ... ){

  # collect data
  xx <- lapply( design$variables[ rounds + 1 ] , function( z ) stats::model.frame( x , data = z , na.action = stats::na.pass ) )
  xx <- do.call( cbind , xx )

  # test column format
  if ( !all( sapply( xx , is.factor ) ) ) stop( "this function is only valid for factors." )
  xlevels <- lapply( xx , function(zz) levels( zz ) )
  ll <- lapply( xlevels , function(zz) length( zz ) )
  if ( !all( ll[[1]] == ll[[2]] && xlevels[[1]] == xlevels[[2]] ) ) stop( "inconsistent categories across rounds." )
  xlevels <- unique( unlist(xlevels) )
  ll <- unique( unlist(ll) )
  has.order <- ifelse( all( sapply( xx , is.ordered ) ) , TRUE , FALSE )

  # add time frame
  colnames( xx ) <- paste0( "round" , seq_along( colnames( xx ) ) - 1 , ":" , colnames( xx ) )

  # collect weights
  ww <- stats::weights( design , "sampling" )

  # aggregate
  NN <- stats::xtabs( c(ww,0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE , drop.unused.levels = FALSE )
  RR <- NN[ , ncol(NN) ][ - nrow( NN ) ]
  CC <- NN[ nrow(NN) , ][ - ncol( NN ) ]
  MM <- NN[ nrow( NN ) , ncol( NN ) ]
  NN <- as.matrix( NN[ -nrow( NN ) , -ncol( NN ) ] )
  N  <- sum( NN ) + sum( RR ) + sum( CC ) + MM

  # handling missing
  if ( sum( MM , RR , CC ) > 0 && !na.rm ) {
    # format results
    NN[,] <- NA
    rval <- NN
    rval <- unclass(rval)
    class(rval) <- "flowstat"
    attr( rval , "var" )       <- unclass(NN)
    attr( rval , "statistic" ) <- flow.type
    attr( rval , "rounds" )    <- rounds
    attr( rval , "formula" )   <- x
    attr( rval , "has.order" )   <- has.order
    return(rval)
  }

  # model fitting
  res <- ipf( xx , ww )

  # collect results
  nipij <- res[["nipij"]]
  eta_i <- res[["eta_i"]]
  p_ij <- res[["p_ij"]]
  rval <- if ( flow.type == "gross" ) N * nipij else nipij
  if ( any( (N * nipij) == 0 ) ) warning( "Some flows are equal to zero. Variances are inconsistent." )

  ### variance calculations

  # get replication weights
  wr <- stats::weights( design , "analysis" )

  # calculate replicates
  lres <- lapply( seq_len(ncol(wr)) , function( irep ) {

    # model fitting
    res_rep <- ipf( xx , wr[ , irep ] )

    # collect res_repults
    nipij_r <- res_rep[["nipij"]]
    eta_i_r <- res_rep[["eta_i"]]
    p_ij_r <- res_rep[["p_ij"]]
    N_r <- res_rep[["N"]]

    list( "gross" = nipij_r * N_r , "net" = nipij_r , "eta_i" = eta_i_r , "p_ij" = p_ij_r )

  } )

  if (extra) {
    # transition probabilities
    qq <- t( Reduce( cbind , lapply( lres , function(zz) matrix( t( zz[["p_ij"]] ) ) ) ) )
    p_var <- diag( survey::svrVar( qq , design$scale , design$rscales , mse = design$mse , coef = matrix( t( res[["p_ij"]] ) ) ) )
    pij_res <- p_ij
    class( pij_res ) <- "flowstat"
    attr( pij_res , "var" ) <- p_var
    attr( pij_res , "statistic" ) <- "transition probabilities"

    # initial probabilities
    qq <- Reduce( rbind , lapply( lres , function(zz) zz[["eta_i"]] ) )
    eta_var <- survey::svrVar( qq , design$scale , design$rscales , mse = design$mse , coef = matrix( t( res[["eta_i"]] ) ) )
    eta_res <- eta_i
    colnames( eta_var ) <- rownames( eta_var ) <-  names( eta_res ) <- xlevels
    class( eta_res ) <- "svystat"
    attr( eta_res , "var" ) <- eta_var
    attr( eta_res , "statistic" ) <- "initial probabilities"
  }

  # decompose list
  reps <- lapply( lres , function(zz) zz[[flow.type]] )
  null_reps <- sapply( reps , function(zz) is.null(zz) )
  reps <- reps[!null_reps]
  qq <- t( Reduce( cbind , lapply( reps , function(zz) matrix( t( zz ) ) ) ) )
  vest <- survey::svrVar( qq , design$scale , design$rscales[ !null_reps ] , mse = design$mse , coef = t( matrix( rval ) ) )
  mu_var <- matrix( diag( vest ) , nrow = nrow( NN ) , ncol = ncol( NN ) , byrow = TRUE )

  # format results
  rval <- unclass(rval)
  class(rval) <- "flowstat"
  attr( rval , "var" ) <- unclass( mu_var )
  dimnames( attr( rval , "var" ) ) <- dimnames( rval )
  attr( rval , "statistic" ) <- flow.type
  attr( rval , "rounds" )    <- rounds
  attr( rval , "formula" )   <- x
  attr( rval , "has.order" )   <- has.order
  attr( rval , "eta" ) <- if (extra) eta_res else NULL
  attr( rval , "pij" ) <- if (extra) pij_res else NULL
  rval

}

#' @export
#' @rdname svyflow
#' @method svyflow surflow.design
svyflow.surflow.design <- function( x , design , flow.type = NULL , rounds = c(0,1) , max.iter = 10 , extra = FALSE , na.rm = FALSE , ... ) {
  flow.type <- match.arg( flow.type , c( "gross" , "net" ) )
  if ( !all( rounds %in% c(0, seq_along( design$variables ) ) ) ) stop( "rounds not in range." )
  NextMethod( "svyflow" , design = design , flow.type = flow.type , rounds = rounds , max.iter = max.iter , extra = extra , na.rm = na.rm , ... )
}

#' @export
svyflow <- function( x , design , ... ) {
  # test valid arguments
  if ( class(x) != "formula" ) stop( "x must be a formula." )
  if ( length(x) != 2 ) stop( "x must be a single variable." )
  UseMethod( "svyflow" , design )
}
