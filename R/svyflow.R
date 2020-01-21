#' @name svyflow
#' @title Summary statistics for repeated sample surveys
#'
#' @description Compute gross flows for data from complex surveys with repeated samples.
#'
#' @param x  A formula
#' @param design  surflow.design object
#' @param flow.type  type of flow to estimate: "gross" for counts and "net" for probabilities. Defaults to \code{flow.type = "gross"}.
#' @param rounds  a vector of integers indicating which round to use. Defaults to \code{rounds = c(0,1)}.
#' @param max.iter  number of iterations. Defaults to \code{max.iter = 10}.
#' @param ...  future expansion.
#'
#' @details ...
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
svyflow.survey.design2 <- function( x , design , flow.type , rounds , max.iter , ... ){

  # preprocess
  design$variables <- lapply( design$variables , function( zz ) {rownames(zz) <- seq_len( nrow(zz)) ; zz} )

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

  # add time frame
  colnames( xx ) <- paste0( "round" , seq_along( colnames( xx ) ) , ":" , colnames( xx ) )

  # collect weights
  ww <- 1 / design$prob

  # aggregate
  NN <- stats::xtabs( c(ww,0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE , drop.unused.levels = FALSE )
  RR <- NN[ , ncol(NN) ][ - nrow( NN ) ]
  CC <- NN[ nrow(NN) , ][ - ncol( NN ) ]
  MM <- NN[ nrow( NN ) , ncol( NN ) ]
  NN <- as.matrix( NN[ -nrow( NN ) , -ncol( NN ) ] )
  N  <- sum( NN ) + sum( RR ) + sum( CC ) + MM

  # maximum pseudo-likelihood estimates for psi, rhoRR, and rhoMM (Rojas et al., 2014, p.296 , Result 4.2 )
  psi <- ( sum( NN ) + sum( RR ) ) / N
  rhoRR <- sum( NN ) / ( sum( NN ) + sum( RR ) )
  rhoMM <- MM / (sum( CC ) + MM )

  ### maximum pseudo-likelihood estimates for eta_i and p_ij (Rojas et al., 2014, p.296 , Result 4.3 )

  # starting values, using Chen and Fienberg (1974) reccomendation
  eta_iv <- rowSums( NN ) / sum( NN )
  p_ijv   <- sweep( NN , 1 , rowSums( NN ) , "/" )

  # iterative process
  v = 0
  while( v < max.iter ) {
    # calculate values
    nipij <- sweep( p_ijv , 1 , eta_iv , FUN = "*" )
    eta_iv <- ( rowSums( NN ) + RR + rowSums( sweep( nipij , 2 , CC / colSums( nipij ) , "*" ) ) ) / ( sum( NN ) + sum( RR ) + sum( CC ) )
    p_ijv <- sweep( NN + sweep( nipij , 2 , CC / colSums( nipij ) , "*" ) , 1:2 , rowSums( NN ) + rowSums( sweep( nipij , 2 , CC / colSums( nipij ) , "*" ) ) , "/" )
    v     <- v + 1
  }
  nipij <- sweep( p_ijv , 1 , eta_iv , FUN = "*" )
  mu_ij <- N * nipij

  ### variance calculations

  # yy array
  yy <- array( 0  , dim = c( nrow( xx ) , ll , ncol( xx ) ) )
  for ( r in seq_len( ncol( xx ) ) ) {
    kk <- stats::model.matrix( ~-1+. , data = xx[,r,drop = FALSE] , contrasts.arg = lapply( xx[,r, drop = FALSE ] , stats::contrasts, contrasts=FALSE ) , na.action = stats::na.pass )
    yy[ match( rownames( kk ) , seq_len( nrow( yy ) ) ) , , r ] <- kk ; rm( kk )
  }

  # create matrix of z variables
  zz <- apply( xx , 2 , function(z) as.numeric( !is.na(z) ) )

  ### special variables
  vv_k <- rowSums( yy[,,1] ) * rowSums( yy[,,2] ) + rowSums( yy[,,2] * (1 - zz[,1]) ) + rowSums( yy[,,1] * (1 - zz[,2]) ) + ( 1- zz[,1] ) * ( 1- zz[,2] )
  nnij_k <- array( 0  , dim = c( nrow( xx ) , nrow( NN ) , ncol( NN ) ) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) nnij_k[,i,j] <- yy[,i,1] * yy[,j,2]

  # u variables
  ueta_ik <- array( 0 , dim = c( nrow(xx) , nrow(NN) ) )
  for ( i in seq_len( nrow(NN) ) ) {
    ueta_ik[,i] <-
      ( rowSums( nnij_k[,i,] ) + yy[,i,1] * ( 1 - zz[,2] ) ) / eta_iv [i] +
      rowSums( sweep( sweep( yy[,,2] , 1 , ( 1 - zz[,1] ) , "*" ) , 2 , p_ijv[i,] / colSums( nipij ) , "*" ) ) +
      (1 - zz[,1]) * (1 - zz[,2])
  }

  up_ijk <- array( NA , dim = c( nrow(xx) , nrow(NN) , ncol(NN) ) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    up_ijk[,i,j] <- nnij_k[,i,j] / p_ijv[i,j] + yy[,i,1] * (1 - zz[,2]) +
      yy[,j,2] * (1 - zz[,1]) * ( eta_iv[i] / colSums( nipij )[j] ) + (1 - zz[,1]) * (1 - zz[,2]) * eta_iv[i]
  }

  # j variables
  jeta_i <- vector( "numeric" , length = nrow(NN) )
  for ( i in seq_len( nrow(NN) ) ) {
    jeta_i[i] <- -2/(eta_iv[i])^2 * sum( ww * yy[,i,1] ) +
      (1/eta_iv[i])^2 * sum( ww * yy[,i,1] * zz[,2] ) - sum( ww * ( 1 - zz[,1] ) * rowSums( sweep( yy[,,2] , 2 , (p_ijv[i,] / colSums( nipij ))^2 , "*" ) ) )
  }

  jp_ij <- array( NA , dim = c( nrow(NN) , ncol(NN) ) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    jp_ij[i,j] <- - 1/(p_ijv[i,j]^2) * sum( ww * nnij_k[,i,j] ) - ( ( eta_iv[i] / colSums( nipij )[j] )^2 ) * sum( ww * yy[,j,2] * (1 - zz[,1] ) )
  }

  # # variance of eta and p-matrix
  # eta_var <- survey::svyrecvar( ww * ueta_ik , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  # eta_var <- diag(eta_var) / jeta_i^2
  # p_var <- matrix( NA, nrow = nrow(NN) , ncol = ncol(NN) )
  # for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
  #   p_var[i,j] <- survey::svyrecvar( ww * up_ijk[[i]][[j]] , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  #   p_var[i,j] <- p_var[i,j] / jp_ij[[i]][[j]]^2
  # }
  #
  # # variance of NN
  # NNvar <- matrix( NA, nrow = nrow(NN) , ncol = ncol(NN) )
  # for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
  #   NNvar[i,j] <- survey::svyrecvar( ww * nnij_k[[i]][[j]] , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  # }
  #
  # # a variables
  # mu_var <- matrix( NA, nrow = nrow(NN) , ncol = ncol(NN) )
  # for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
  #   a7 <- nipij[i,j]
  #   a8 <- NN[i,j] * p_ijv[i,j]
  #   a9 <- NN[i,j] * eta_iv[i]
  #   mu_var[i,j] <- (a7^2) * NNvar[i,j] + (a8^2) * eta_var[i] + (a9^2) * p_var[i,j]
  # }

  # alternative strategy
  nipij_k <- array( NA, dim = c( nrow( xx ) , nrow( NN ) , ncol( NN ) ) )
  muij_k <- array( NA, dim = c( nrow( xx ) , nrow( NN ) , ncol( NN ) ) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    nipij_k[,i,j] <- p_ijv[i,j] * ( ueta_ik[,i] / jeta_i[i] ) + eta_iv[i] * ( up_ijk[,i,j] / jp_ij[i,j] )
    muij_k[,i,j] <- nipij[i,j] * vv_k + N * nipij_k[,i,j]
  }
  mu_var <- matrix( NA, nrow = nrow(NN) , ncol = ncol(NN) )
  for ( i in seq_len( nrow(NN) ) ) {
    if ( flow.type == "gross" ) {
      mu_var[i,] <- diag( survey::svyrecvar( ww * muij_k[,i,] , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    } else {
      mu_var[i,] <- diag( survey::svyrecvar( ww * nipij_k[,i,] , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    }
  }

  # format results
  rval <- if ( flow.type == "gross" ) mu_ij else nipij
  rval <- as.table(rval)
  class(rval) <- "flowstat"
  attr( rval , "var" )       <- as.table(mu_var)
  dimnames( attr( rval , "var" ) ) <- dimnames( mu_ij )
  attr( rval , "statistic" ) <- flow.type
  attr( rval , "rounds" )    <- rounds
  attr( rval , "formula" )   <- x
  rval

}

#' @export
#' @rdname svyflow
#' @method svyflow svyrep.design
svyflow.svyrep.design <- function( x , design , flow.type , rounds , max.iter , ... ){

  # collect data
  xx <- lapply( design$variables[ rounds + 1 ] , function( z ) stats::model.frame( x , data = z , na.action = stats::na.pass ) )
  xx <- do.call( cbind , xx )

  # test column format
  if ( !all( sapply( xx , is.factor ) ) ) stop( "this function is only valid for factors." )

  # add time frame
  colnames( xx ) <- paste0( "round" , seq_along( colnames( xx ) ) , ":" , colnames( xx ) )

  # collect weights
  ww <- stats::weights( design , "sampling" )

  # aggregate
  NN <- stats::xtabs( c(ww,0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE , drop.unused.levels = FALSE )
  RR <- NN[ , ncol(NN) ][ - nrow( NN ) ]
  CC <- NN[ nrow(NN) , ][ - ncol( NN ) ]
  MM <- NN[ nrow( NN ) , ncol( NN ) ]
  NN <- as.matrix( NN[ -nrow( NN ) , -ncol( NN ) ] )
  N  <- sum( NN ) + sum( RR ) + sum( CC ) + MM

  # maximum pseudo-likelihood estimates for psi, rhoRR, and rhoMM (Rojas et al., 2014, p.296 , Result 4.2 )
  psi <- ( sum( NN ) + sum( RR ) ) / N
  rhoRR <- sum( NN ) / ( sum( NN ) + sum( RR ) )
  rhoMM <- MM / (sum( CC ) + MM )

  ### maximum pseudo-likelihood estimates for eta_i and p_ij (Rojas et al., 2014, p.296 , Result 4.3 )

  # starting values, using Chen and Fienberg (1974) reccomendation
  eta_iv <- rowSums( NN ) / sum( NN )
  p_ijv   <- sweep( NN , 1 , rowSums( NN ) , "/" )

  # iterative process
  v = 0
  while( v < max.iter ) {
    # calculate values
    nipij <- sweep( p_ijv , 1 , eta_iv , FUN = "*" )
    eta_iv <- ( rowSums( NN ) + RR + rowSums( sweep( nipij , 2 , CC / colSums( nipij ) , "*" ) ) ) / ( sum( NN ) + sum( RR ) + sum( CC ) )
    p_ijv <- sweep( NN + sweep( nipij , 2 , CC / colSums( nipij ) , "*" ) , 1:2 , rowSums( NN ) + rowSums( sweep( nipij , 2 , CC / colSums( nipij ) , "*" ) ) , "/" )
    v     <- v + 1
  }
  nipij <- sweep( p_ijv , 1 , eta_iv , FUN = "*" )
  rval <- if ( flow.type == "gross" ) N * nipij else nipij

  ### variance calculations

  # get replication weights
  wr <- stats::weights( design , "analysis" )

  # calculate replicates
  lres <- lapply( seq_len(ncol(wr)) , function( zz ) {

    # aggregate
    NN <- stats::xtabs( c(wr[,zz],0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE , drop.unused.levels = FALSE )
    RR <- NN[ , ncol(NN) ][ - nrow( NN ) ]
    CC <- NN[ nrow(NN) , ][ - ncol( NN ) ]
    MM <- NN[ nrow( NN ) , ncol( NN ) ]
    NN <- as.matrix( NN[ -nrow( NN ) , -ncol( NN ) ] )
    N  <- sum( NN ) + sum( RR ) + sum( CC ) + MM

    # maximum pseudo-likelihood estimates for psi, rhoRR, and rhoMM (Rojas et al., 2014, p.296 , Result 4.2 )
    psi <- ( sum( NN ) + sum( RR ) ) / N
    rhoRR <- sum( NN ) / ( sum( NN ) + sum( RR ) )
    rhoMM <- MM / (sum( CC ) + MM )

    ### maximum pseudo-likelihood estimates for eta_i and p_ij (Rojas et al., 2014, p.296 , Result 4.3 )

    # starting values, using Chen and Fienberg (1974) reccomendation
    eta_iv <- rowSums( NN ) / sum( NN )
    p_ijv   <- sweep( NN , 1 , rowSums( NN ) , "/" )

    # iterative process
    v = 0
    while( v < max.iter ) {
      # calculate values
      nipij <- sweep( p_ijv , 1 , eta_iv , FUN = "*" )
      eta_iv <- ( rowSums( NN ) + RR + rowSums( sweep( nipij , 2 , CC / colSums( nipij ) , "*" ) ) ) / ( sum( NN ) + sum( RR ) + sum( CC ) )
      p_ijv <- sweep( NN + sweep( nipij , 2 , CC / colSums( nipij ) , "*" ) , 1:2 , rowSums( NN ) + rowSums( sweep( nipij , 2 , CC / colSums( nipij ) , "*" ) ) , "/" )
      v     <- v + 1
    }
    nipij <- sweep( p_ijv , 1 , eta_iv , FUN = "*" )
    mu_ij <- N * nipij

    list( "gross" = mu_ij , "net" = nipij )

  } )

  # decompose list
  reps <- lapply( lres , function(zz) zz[[flow.type]] )
  null_reps <- sapply( reps , function(zz) is.null(zz) )
  reps <- reps[!null_reps]
  p_var <- NN
  p_var[,] <- NA
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    qq <- sapply( reps , function(zz) zz[ i , j ] )
    p_var[i,j] <- survey::svrVar( qq , design$scale , design$rscales[ !null_reps ] , mse = design$mse , coef = rval[i,j] )
  }

  # format results
  rval <- as.table(rval)
  class(rval) <- "flowstat"
  attr( rval , "var" ) <- p_var
  dimnames( attr( rval , "var" ) ) <- dimnames( rval )
  attr( rval , "statistic" ) <- flow.type
  attr( rval , "rounds" )    <- rounds
  attr( rval , "formula" )   <- x
  rval

}

#' @export
#' @rdname svyflow
#' @method svyflow surflow.design
svyflow.surflow.design <- function( x , design , flow.type = NULL , rounds = c(0,1) , max.iter = 10 , ... ) {
  flow.type <- match.arg( flow.type , c( "gross" , "net" ) )
  if ( !all( rounds %in% c(0, seq_along( design$variables ) ) ) ) stop( "rounds not in range." )
  NextMethod( "svyflow" , design = design , flow.type = flow.type , rounds = rounds , max.iter = max.iter , ... )
}

#' @export
svyflow <- function( x , design , ... ) {
  # test valid arguments
  if ( class(x) != "formula" ) stop( "x must be a formula." )
  if ( length(x) != 2 ) stop( "x must be a single variable." )
  UseMethod( "svyflow" , design )
}
