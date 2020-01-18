#' @name svyflow
#' @title Summary statistics for repeated sample surveys
#'
#' @description Compute gross flows for data from complex surveys with repeated samples.
#'
#' @param x  A formula
#' @param design  survflow.design object
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
#' @references Gutierrez, A., Trujillo, L. \& Silva, N. (2014). The estimation of gross flows in complex surveys with random nonresponse,
#' Survey Methodology 40(2), pp. 285-321.
#'
#' Lumley, Thomas S. (2010). Complex Surveys: A Guide to Analysis Using R. Wiley Publishing.
#'
#' @examples
#' # load data
#' data( "artificial" )
#'
#' # create surf design object
#' flowdes <-
#'   svyflowdesign( ids = ~0 ,
#'                  probs = ~ prob ,
#'                  data.list = list( dfa0 , dfa1 ) ,
#'                  nest = TRUE )
#'
#' # gross flows
#' gross.flows <- svyflow( ~v0 , design = flowdes , flow.type = "gross" )
#' SE( gross.flows )
#'
#' # net flows
#' net.flows <- svyflow( ~v0 , design = flowdes , flow.type = "net" )
#' SE( net.flows )
#'
#' @export
#' @rdname svyflow
#' @method svyflow survey.design2
svyflow.survey.design2 <- function( x , design , flow.type , rounds , max.iter , ... ){

  # collect data
  xx <- lapply( design$variables[ rounds + 1 ] , function( z ) stats::model.frame( x , data = z , na.action = na.pass ) )
  xx <- do.call( cbind , xx )

  # test column format
  if ( !all( sapply( xx , is.factor ) ) ) stop( "this function is only valid for factors." )

  # add time frame
  colnames( xx ) <- paste0( "round" , seq_along( colnames( xx ) ) , ":" , colnames( xx ) )

  # collect weights
  ww <- 1 / design$prob

  # aggregate
  NN <- xtabs( c(ww,0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE )
  RR <- NN[ , ncol(NN) ][ - nrow( NN ) ]
  CC <- NN[ nrow(NN) , ][ - ncol( NN ) ]
  MM <- NN[ nrow( NN ) , ncol( NN ) ]
  NN <- as.matrix( NN[ -nrow( NN ) , -ncol( NN ) ] )
  N  <- sum( NN ) + sum( RR ) + sum( CC ) + MM
  if ( any( NN <= 0 ) ) stop( "Some flows are equal to zero.")

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

  # yy matrix
  yy <- lapply( seq_len( ncol( xx ) ) , function( j ) {
    kk <- model.matrix( ~-1+. , data = xx[,j,drop = FALSE] , contrasts.arg = lapply( xx[ , j, drop = FALSE ] , contrasts, contrasts=FALSE ) , na.action = na.pass )
    oo <- matrix( 0 , nrow = nrow(xx) , ncol = ncol(kk) , dimnames = list( seq_len( nrow(xx) ) , colnames( kk) ) )
    oo[ match( rownames( kk ) , rownames( oo )) , ] <- kk
    oo
  } )

  # create matrix of contrasts and z variables
  zz <- apply( xx , 2 , function(z) as.numeric( !is.na(z) ) )

  ### special variables
  vv_k <- rowSums( yy[[1]] ) * rowSums( yy[[2]] ) + rowSums( yy[[2]] * (1 - zz[,1]) ) + rowSums( yy[[1]] * (1 - zz[,2]) )
  nnij_k <- list(NULL)
  for ( i in seq_len( nrow(NN) ) ) nnij_k[[i]] <- list(NULL)
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) nnij_k[[i]][[j]] <- yy[[1]][,i] * yy[[2]][,j]

  # u variables
  ueta_ik <- list(NULL)
  for ( i in seq_len( nrow(NN) ) ) {
    ueta_ik[[i]] <-
      ( rowSums( do.call( cbind , nnij_k[[i]] ) ) + yy[[1]][,i, drop = FALSE ] * ( 1 - zz[,2 , drop = FALSE ] ) ) / eta_iv [i] +
      rowSums( sweep( sweep( yy[[2]] , 1 , ( 1 - zz[,1] ) , "*" ) , 2 , p_ijv[i,] / colSums( nipij ) , "*" ) ) +
      (1 - zz[,1]) * (1 - zz[,2])
  }
  ueta_ik <- do.call( cbind , ueta_ik)

  up_ijk <- list(NULL)
  for ( i in seq_len( nrow(NN) ) ) up_ijk[[i]] <- list(NULL)
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    up_ijk[[i]][[j]] <- nnij_k[[i]][[j]] / p_ijv[i,j] + yy[[1]][,i] * (1 - zz[,2]) +
      yy[[2]][,j] * (1 - zz[,1]) * ( eta_iv[i] / colSums( nipij )[j] ) + (1 - zz[,1]) * (1 - zz[,2]) * eta_iv[i]
  }

  # j variables
  jeta_i <- list(NULL)
  for ( i in seq_len( nrow(NN) ) ) {
    jeta_i[[i]] <- -2/(eta_iv[i])^2 * sum( ww * yy[[1]][,i] ) +
      (1/eta_iv[i])^2 * sum( ww * yy[[1]][,i] * zz[,2] ) - sum( ww * ( 1 - zz[,1] ) * rowSums( sweep( yy[[2]] , 2 , (p_ijv[i,] / colSums( nipij ))^2 , "*" ) ) )
  }
  jeta_i <- do.call( cbind , jeta_i )

  jp_ij <- list(NULL)
  for ( i in seq_len( nrow(NN) ) ) jp_ij[[i]] <- list(NULL)
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    jp_ij[[i]][[j]] <- - 1/(p_ijv[i,j]^2) * sum( ww * nnij_k[[i]][[j]] ) - ( ( eta_iv[i] / colSums( nipij )[j] )^2 ) * sum( ww * yy[[2]][,j] * (1 - zz[,1] ) )
  }

  # variance of NN
  NNvar <- matrix( NA, nrow = nrow(NN) , ncol = ncol(NN) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    NNvar[i,j] <- survey::svyrecvar( ww * nnij_k[[i]][[j]] , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  }

  # variance of eta and p-matrix
  eta_var <- survey::svyrecvar( ww * ueta_ik , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  eta_var <- diag(eta_var) / jeta_i^2
  p_var <- matrix( NA, nrow = nrow(NN) , ncol = ncol(NN) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    p_var[i,j] <- survey::svyrecvar( ww * up_ijk[[i]][[j]] , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    p_var[i,j] <- p_var[i,j] / jp_ij[[i]][[j]]^2
  }

  # a variables
  mu_var <- matrix( NA, nrow = nrow(NN) , ncol = ncol(NN) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    a7 <- nipij[i,j]
    a8 <- NN[i,j] * p_ijv[i,j]
    a9 <- NN[i,j] * eta_iv[i]
    mu_var[i,j] <- (a7^2) * NNvar[i,j] + (a8^2) * eta_var[i] + (a9^2) * p_var[i,j]
  }

  # format results
  rval <- if ( flow.type == "gross" ) mu_ij else nipij
  class(rval) <- "flowstat"
  attr( rval , "var" )       <- if ( flow.type == "gross" ) mu_var else mu_var / NN^2
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
  stop( "replicate designs are not supported." )
}

#' @export
#' @rdname svyflow
#' @method svyflow survflow.design
svyflow.survflow.design <- function( x , design , flow.type = c( "gross" , "net" ) , rounds = c(0,1) , max.iter = 10 , ... ) {
  flow.type <- match.arg( flow.type )
  if ( !all( rounds %in% c(0, seq_along( design$variables ) ) ) ) stop( "rounds not in range." )
  NextMethod( "svyflow" , design = design , flowtype , rounds , max.iter , ... )
}

#' @export
svyflow <- function( x , design , ... ) {
  # test valid arguments
  if ( class(x) != "formula" ) stop( "x must be a formula." )
  if ( length(x) != 2 ) stop( "x must be a single variable." )
  UseMethod( "svyflow" , design )
}
