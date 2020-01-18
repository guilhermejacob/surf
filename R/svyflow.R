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
#' These objects have methods for coef , vcov, SE, and cv.
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
#' data( "initial" )
#' data( "final" )
#'
#' # create survflow design object
#' flowdes <-
#'   svyflowdesign( ids = ~ upa ,
#'                  strata = ~ estrato ,
#'                  probs = ~ longprob ,
#'                  data.list = list( f.qtr , s.qtr ) ,
#'                  nest = TRUE )
#'
#' # totals
#' svytotal( ~factor( vd4002 ) , flowdes , na.rm = TRUE ) # drops observations with missing in any rounds
#' svytotal( ~factor( vd4002 ) , flowdes , na.rm = TRUE , which.dataset = 0 ) # drops observations with missing in the first round
#' svytotal( ~factor( vd4002 ) , flowdes , na.rm = TRUE , which.dataset = 1 ) # drops observations with missing in the second round
#'
#' # means
#' svymean( ~factor( vd4002 ) , flowdes , na.rm = TRUE ) # drops observations with missing in any rounds
#' svymean( ~factor( vd4002 ) , flowdes , na.rm = TRUE , which.dataset = 0 ) # drops observations with missing in the first round
#' svymean( ~factor( vd4002 ) , flowdes , na.rm = TRUE , which.dataset = 1 ) # drops observations with missing in the second round
#'
#' @export
#' @rdname svyflow
#' @method svyflow survey.design2
svyflow.survey.design2 <- function( x , design , flow.type = c( "gross" , "net" ) , rounds = c(0,1) , max.iter = 10 , ... ){

  # test valid arguments
  flow.type <- match.arg( flow.type )
  if ( class(x) != "formula" ) stop( "x must be a formula." )
  if ( length(rounds) != 2 ) stop( "rounds must have length = 2." )
  if ( !all( rounds %in% c(0, seq_along( design$variables ) ) ) ) stop( "rounds not in range." )

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
  NN <- xtabs( ww ~ . , data = xx , addNA = TRUE )
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
  xxf <- model.matrix( ~-1+. , data = xx , contrasts.arg = lapply( xx , contrasts, contrasts=FALSE ) , na.action = na.pass )
  oo <- matrix( 0 , nrow = nrow(xx) , ncol = ncol(xxf) , dimnames = list( seq_len( nrow(xx) ) , NULL ) )
  oo[ match( rownames( xxf ) , rownames( oo ) ) , ] <- xxf
  xxf <- oo ; rm( oo )
  zz <- apply( xx , 2 , function(z) as.numeric( !is.na(z) ) )

  ### special variables
  vv_k <- rowSums( yy[[1]] ) * rowSums( yy[[2]] ) + rowSums( yy[[2]] * (1 - zz[,1]) ) + rowSums( yy[[1]] * (1 - zz[,2]) )
  nnij_k <- lapply( seq_len( ncol( yy[[1]] ) ) , function(z) yy[[1]][ , z ] * yy[[2]] )

  # u variables
  ueta_ik <- list(NULL)
  for ( i in seq_len( nrow(NN) ) ) {
    ueta_ik[[i]] <-
      ( rowSums( nnij_k[[i]] ) + yy[[1]][,i, drop = FALSE ] * ( 1 - zz[,2 , drop = FALSE ] ) ) / eta_iv [ i ] +
      rowSums( sweep( yy[[2]] * ( 1 - zz[,1] ) , 2 , p_ijv[i,] / colSums( nipij ) , "*" ) ) +
      (1 - zz[,1]) * (1 - zz[,2])
  }
  ueta_ik <- do.call( cbind , ueta_ik)

  up_ijk <- list(NULL)
  for ( i in seq_len( nrow(NN) ) ) up_ijk[[i]] <- list(NULL)
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    up_ijk[[i]][[j]] <- nnij_k[[i]][,j] / p_ijv[i,j] + yy[[1]][,i] * (1 - zz[,2]) +
      yy[[2]][,j] * (1 - zz[,1]) * ( eta_iv[i] * colSums( nipij )[j] ) + (1 - zz[,1]) * (1 - zz[,2]) * eta_iv[i]
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
    jp_ij[[i]][[j]] <- - (1/p_ijv[i,j])^2 * sum( ww * yy[[1]][ , i ] * yy[[2]][ , j ] ) - ( eta_iv[i] / rowSums( nipij^2 )[j] ) * sum( ww * yy[[2]][,j] * (1 - zz[,1]) )
  }

  # corrects u variables; i.e., u / j
  ueta_ik <- sweep( ueta_ik , 2 , jeta_i ,  "/" )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    up_ijk[[i]][[j]] <- up_ijk[[i]][[j]] / jp_ij[[i]][[j]]
  }

  # variance of NN
  NNvar <- matrix( NA, nrow = nrow(NN) , ncol = ncol(NN) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    NNvar[i,j] <- survey::svyrecvar( nnij_k[[i]][,j] * ww , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  }

  # variance of eta and p-matrix
  eta_var <- survey::svyrecvar( ueta_ik * ww , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  p_var <- matrix( NA, nrow = nrow(NN) , ncol = ncol(NN) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    p_var[i,j] <- survey::svyrecvar( up_ijk[[i]][[j]] * ww , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  }

  # a variables
  mu_var <- matrix( NA, nrow = nrow(NN) , ncol = ncol(NN) )
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    a7 <- nipij[i,j]
    a8 <- nipij[i,j]
    a9 <- nipij[i,j]
    mu_var[i,j] <- a7^2 * NNvar[i,j] + a8^2 * eta_var[i,i] + a9^2 * p_var[i,j]
  }

  # format results
  rval <- if ( flow.type == "gross" ) mu_ij else nipij
  class(rval) <- "flowstat"
  attr( rval , "var" )       <- if ( flow.type == "gross" ) mu_var else mu_var / N^2
  dimnames( attr( rval , "var" ) ) <- dimnames( mu_ij )
  attr( rval , "statistic" ) <- flow.type
  attr( rval , "rounds" )    <- rounds
  attr( rval , "formula" )   <- x
  rval

}

#' @export
#' @rdname svyflow
#' @method svyflow svyrep.design
svyflow.svyrep.design <- function( x , design , flow.type = c( "gross" , "net" ) , rounds = c(0,1) , max.iter = 10 , ... ){

  # test valid arguments
  flow.type <- match.arg( flow.type )
  if ( class(x) != "formula" ) stop( "x must be a formula." )
  if ( length(rounds) != 2 ) stop( "rounds must have length = 2." )
  if ( !all( rounds %in% c(0, seq_along( design$variables ) ) ) ) stop( "rounds not in range." )

  # collect data
  xx <- lapply( design$variables[ rounds + 1 ] , function( z ) stats::model.frame( x , data = z , na.action = na.pass ) )
  xx <- do.call( cbind , xx )

  # test column format
  if ( !all( sapply( xx , is.factor ) ) ) stop( "this function is only valid for factors." )

  # add time frame
  colnames( xx ) <- paste0( "round" , seq_along( colnames( xx ) ) , ":" , colnames( xx ) )

  # collect weights
  ww <- weights( design , "sampling" )

  # aggregate
  NN <- xtabs( ww ~ . , data = xx , addNA = TRUE )
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
  rval <- if ( flow.type == "gross" ) N * nipij else nipij

  ### variance calculations

  # get replication weights
  wr <- weights( design , "analysis" )

  # calculate replicates
  lres <- lapply( seq_len(ncol(wr)) , function( zz ) {

    # aggregate
    NN <- xtabs( wr[,zz] ~ . , data = xx , addNA = TRUE )
    RR <- NN[ , ncol(NN) ][ - nrow( NN ) ]
    CC <- NN[ nrow(NN) , ][ - ncol( NN ) ]
    MM <- NN[ nrow( NN ) , ncol( NN ) ]
    NN <- as.matrix( NN[ -nrow( NN ) , -ncol( NN ) ] )
    N  <- sum( NN ) + sum( RR ) + sum( CC ) + MM
    if ( any( NN <= 0 ) ) return( NULL )

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
  p_var <- rval
  p_var[,] <- NA
  for ( i in seq_len( nrow(NN) ) ) for ( j in seq_len( ncol( NN ) ) ) {
    qq <- sapply( reps , function(zz) zz[ i , j] )
    p_var[i,j] <- survey::svrVar( qq , design$scale , design$rscales[ !null_reps ] , mse = design$mse , coef = rval )
  }

  # format results
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
#' @method svyflow survflow.design
svyflow.survflow.design <- function( x , design , ... ) NextMethod( "svyflow" , design )

#' @export
svyflow <- function( x , design , ... ) UseMethod( "svyflow" , design )
