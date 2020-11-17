# fit flows with complete response
frf <- function( Nij ) {

  # remove non-response counts
  N <- sum( Nij )

  # calculate eta
  eta <- rowSums( Nij ) / N

  # calculate pij
  pij <- sweep( Nij , 1 , rowSums( Nij ) , "/" )

  # calculate muij
  muij <- N * sweep( pij , 1 , eta , "*" )

  # add missing to non-response mechanism
  psi <- NA
  rho <- NA
  tau <- NA

  # build result
  mfit <-
    list(
      "model" = "Full Response" ,
      "iter" = NA ,
      "N" = N ,
      "Nij" = Nij ,
      "psi" = NA ,
      "rho" = NA ,
      "tau" = NA ,
      "eta" = eta ,
      "pij" = pij ,
      "muij" = N * sweep( pij , 1 , eta , "*" ) ,
      "gamma" = colSums( sweep( pij , 1 , eta , "*" ) ) )

  # fill model-related estimates
  mfit[["ll"]] <- NA
  mfit[["observed.counts"]] <- NA
  mfit[["estimated.counts"]] <- NA
  mfit[["observed.props"]] <- NA
  mfit[["chimat"]] <- NA

  # return fit
  return( mfit )

}

# calculate variance flows with complete response
frf_variance <- function( xx , ww , res , design ) {

  # load objects
  Nij <- res[["Nij"]]
  N <- sum(Nij)
  eta <- res[["eta"]]
  pij <- res[["pij"]]
  muij <- res[["muij"]]
  nipij <- sweep( pij , 1 , eta , "*" )

  # linearization of sum_i Nij

  # yy array - see Rojas et al. (2014, p.294)
  yy <- array( 0  , dim = c( nrow( xx ) , nrow( Nij ) , ncol( xx ) ) )
  for ( r in seq_len( ncol( xx ) ) ) {
    kk <- stats::model.matrix( ~-1+. , data = xx[,r,drop = FALSE] , contrasts.arg = lapply( xx[,r, drop = FALSE ] , stats::contrasts, contrasts=FALSE ) , na.action = stats::na.pass )
    yy[ which( !is.na( xx[ , r ] ) ) , , r ] <- kk ; rm( kk )
  }

  # Special variables - see Rojas et al. (2014, p.295)
  y1y2 <- array( 0  , dim = c( nrow( xx ) , nrow( Nij ) , ncol( Nij ) ) )
  for ( i in seq_len( nrow( Nij ) ) ) for ( j in seq_len( ncol( Nij ) ) ) y1y2[,i,j] <- yy[,i,1] * yy[,j,2]

  # linearized sum_j Nij
  y1 <- apply( y1y2 , c(1,2) , sum )

  # calculate linearized eta
  u_eta <- sweep( N * y1 , 2 , rowSums( Nij ) , "-" ) / N^2

  # calculate linearized pij
  u_pij <- y1y2
  for ( i in seq_len( nrow( Nij ) ) ) for ( j in seq_len( ncol( Nij ) ) ) u_pij[,i,j] <- ( y1y2[,i,j] * sum( Nij[,j] ) - Nij[i,j] ) / sum( Nij[,j] )^2

  # calculate linearized muij
  u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( Nij ) , ncol( Nij ) ) )
  for ( i in seq_len( nrow(Nij) ) ) for ( j in seq_len( ncol( Nij ) ) ) {
    u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
  }

  # calculate linearized gamma
  u_gamma <- sapply( seq_len( ncol( pij ) ) , function(j) {
    rowSums( sweep( u_pij[,,j] , 2 , eta , "*" ) + sweep( u_eta , 2 , pij[,j] , "*" ) )
  } )

  ### estimate variances

  # create matrix of linearized variables
  u_muij <- do.call( cbind , lapply( seq_len(ncol( Nij)) , function( z ) u_muij[,z,] ) )
  u_pij <- do.call( cbind , lapply( seq_len(ncol( Nij)) , function( z ) u_pij[,z,] ) )

  # fill non-response objects
  psi_var <- as.matrix( NA , length(psi))
  rho_var <- as.matrix( NA , length(rho))
  tau_var <- as.matrix( NA , length(tau))

  # calculate variance of eta
  eta_var <- survey::svyrecvar( sweep( u_eta , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

  # calculate variance of pij
  pij_var <- survey::svyrecvar( ww * u_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  pij_var <- diag( pij_var )
  pij_var <- matrix( pij_var , nrow = nrow( Nij ) , ncol = ncol( Nij ) , byrow = TRUE )

  # calculate variance of muij
  muij_var <- survey::svyrecvar( ww * u_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  muij_var <- diag( muij_var )
  muij_var <- matrix( muij_var , nrow = nrow( Nij ) , ncol = ncol( Nij ) , byrow = TRUE )

  # calculate variance of gamma
  gamma_var <- survey::svyrecvar( sweep( u_gamma , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  # gamma_var <- diag( gamma_var )

  # variance of delta
  u_delta <- sweep( N * ( u_gamma - u_eta ) , 2 , res[["gamma"]] - res[["eta"]] , "-" )
  delta_var <- survey::svyrecvar( sweep( u_delta , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

  ### combine variance estimates
  mvar <-
    list(
      "psi" = psi_var ,
      "rho" = rho_var ,
      "tau" = tau_var ,
      "eta" = eta_var ,
      "pij" = pij_var ,
      "muij" = muij_var ,
      "gamma" = gamma_var ,
      "delta" = delta_var )

  # return object
  return( mvar )

}
