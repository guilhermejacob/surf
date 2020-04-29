# function for model variance calculation
ipf_variance <- function( xx , ww , res , design ) {

  # load objects
  for ( this_obj in names( res ) ) assign( this_obj , res[[this_obj]] )

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

  # ajusts using onde of the models
  mvar <- switch( res$model , MCAR = {

    ### psi

    psi_var <- NA

    ### rhoRR

    rhoRR_var <- NA

    ### rhoMM

    rhoMM_var <- NA

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

    # create matrix of linearized variables
    lin_pij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_pij[,z,] ) )

    # adjust for structural zeros in transition matrix
    if ( !is.null( pij.zero ) ) {
      lin_pij[ , c( t( pij.zero ) ) ] <- 0
    }

    # calculate variance of u_pij sum
    pij_var <- diag( survey::svyrecvar( ww * lin_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    pij_var <- matrix( pij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    # create matrix of linearized variables
    lin_muij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_muij[,z,] ) )

    # adjust for structural zeros in transition matrix
    if ( !is.null( pij.zero ) ) {
      lin_muij[ , c( t( pij.zero ) ) ] <- 0
    }

    # calculate variance of muij
    muij_var <- diag( survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    muij_var <- matrix( muij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### combine results

    # return final estimates
    resvar <-
      list(
        "psi" = psi_var ,
        "rhoRR" = rhoRR_var ,
        "rhoMM" = rhoMM_var ,
        "eta" = eta_var ,
        "pij" = pij_var ,
        "muij" = muij_var )

    # return list of results
    return( resvar )

  } , A = {

    # calculate auxiliary stats
    nipij <- sweep( pij , 1 , eta , "*" )

    ### psi

    # special variables
    a1 <- ( sum( bigCj ) + bigM ) / N^2
    a2 <- ( sum( bigNij ) + sum( bigRi ) ) / N^2

    # linearized variable e_psi
    e_psi <- a1 * ( 2 - zz[,2] ) + a2 * ( 1 - zz[,1] ) * ( 2 - zz[,2] )

    # calculate variance
    psi_var <- survey::svyrecvar( ww * e_psi , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### rhoRR

    # special variables
    a3 <- sum( bigRi ) / ( sum( bigNij ) + sum( bigRi ) )^2
    a4 <- - sum( bigNij ) / ( sum( bigNij ) + sum( bigRi ) )^2

    # linearized variable e_psi
    e_rhoRR <- a3 + a4 * ( 1 - zz[,2] )

    # calculate variance
    rhoRR_var <- survey::svyrecvar( ww * e_rhoRR , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### rhoMM

    # special variables
    a5 <- - bigM / ( sum( bigCj ) + bigM )^2
    a6 <- - sum( bigCj ) / ( sum( bigCj ) + bigM )^2

    # linearized variable e_psi
    e_rhoMM <- a5 * ( 1 - zz[,1] ) + a6 * ( 1 - zz[,1] ) * ( 1 - zz[,2] )

    # calculate variance
    rhoMM_var <- survey::svyrecvar( ww * e_rhoMM , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### eta

    # Calculate scores for estimating the variance of eta parameters
    u_eta <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_eta[,i] <-
        ( rowSums( y1y2[,i,] ) + yy[,i,1] * ( 1 - zz[,2] ) ) / eta[i] +
        rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , pij[i,] / colSums( nipij ) , "*" ) ) +
        ( 1- zz[,1] ) * ( 1 - zz[,2] )
    }

    # Calculate jacobian for estimating the variance of eta parameters
    jeta <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jeta[i] <-
        (-2 / eta[i]^2 ) * sum( ww * yy[,i,1] ) +
        ( 1 / eta[i]^2 ) * sum( ww * yy[,i,1] * zz[,2] ) -
        sum( ww * ( 1 - zz[,1] ) * rowSums( sweep( yy[,,2] , 2 , ( pij[i,] / colSums( nipij ) )^2 , "*" ) ) )
    }

    # divide u_eta by the jacobian
    u_eta <- sweep( u_eta , 2 , jeta , "/" )

    # calculate variance of eta
    eta_var <- survey::svyrecvar( sweep( u_eta , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    eta_var <- diag( eta_var )

    ### pij

    # Calculate scores for estimating the variance of pij parameters
    u_pij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_pij[,i,j] <- ( y1y2[,i,j] / pij[i,j] ) +
        ( yy[,i,1] * ( 1 - zz[,2] ) ) +
        ( yy[,j,2] * ( 1 - zz[,1] ) ) * ( eta[i] / colSums( nipij )[j] ) +
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * eta[i]
    }

    # Calculate jacobian for estimating the variance of pij parameters
    jpij <- array( 0 , dim = c( nrow(bigNij) , ncol(bigNij) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      jpij[i,j] <- ( -1 / pij[i,j]^2 ) * sum( ww * y1y2[,i,j] ) -
        ( ( eta[i] / colSums( nipij )[j] )^2 ) *
        sum( ww * yy[,j,2] * ( 1 - zz[,1] ) )
    }

    # divide u_pij by the jacobian
    u_pij <- sweep( u_pij , 2:3 , jpij , "/" )

    # create matrix of linearized variables
    lin_pij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_pij[,z,] ) )

    # adjust for structural zeros in transition matrix
    if ( !is.null( pij.zero ) ) {
      lin_pij[ , c( t( pij.zero ) ) ] <- 0
    }

    # calculate variance of u_pij sum
    pij_var <- diag( survey::svyrecvar( ww * lin_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    pij_var <- matrix( pij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    # create matrix of linearized variables
    lin_muij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_muij[,z,] ) )

    # adjust for structural zeros in transition matrix
    if ( !is.null( pij.zero ) ) {
      lin_muij[ , c( t( pij.zero ) ) ] <- 0
    }

    # calculate variance of muij
    muij_var <- diag( survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    muij_var <- matrix( muij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### combine results

    # return final estimates
    resvar <-
      list(
        "psi" = psi_var ,
        "rhoRR" = rhoRR_var ,
        "rhoMM" = rhoMM_var ,
        "eta" = eta_var ,
        "pij" = pij_var ,
        "muij" = muij_var )

    # return list of results
    return( resvar )

  } , B ={
    # calculate auxiliary stats
    nipij <- sweep( pij , 1 , eta , "*" )
    psicni <- ( 1 - psi ) * eta
    psicnipij <- sweep( pij , 1 , psicni , "*" )
    psicpij <- sweep( pij , 1 , 1 - psi , "*" )

    ### rhoRR

    # special variables
    a3 <- sum( bigRi ) / ( sum( bigNij ) + sum( bigRi ) )^2
    a4 <- - sum( bigNij ) / ( sum( bigNij ) + sum( bigRi ) )^2

    # linearized variable e_psi
    e_rhoRR <- a3 + a4 * ( 1 - zz[,2] )

    # calculate variance
    rhoRR_var <- survey::svyrecvar( ww * e_rhoRR , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### rhoMM

    # special variables
    a5 <- - bigM / ( sum( bigCj ) + bigM )^2
    a6 <- - sum( bigCj ) / ( sum( bigCj ) + bigM )^2

    # linearized variable e_psi
    e_rhoMM <- a5 * ( 1 - zz[,1] ) + a6 * ( 1 - zz[,1] ) * ( 1 - zz[,2] )

    # calculate variance
    rhoMM_var <- survey::svyrecvar( ww * e_rhoMM , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### psi

    # Calculate scores for estimating the variance of psi parameters
    u_psi <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_psi[,i] <-
        rowSums( y1y2[,i,] ) / psi[i] +
        yy[,i,1] * ( 1 - zz[,2] ) / psi[i] -
        rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , nipij[i,] / colSums( psicnipij ) , "*" ) ) -
        # ( 1- zz[,1] ) * ( 1 - zz[,2] ) * rowSums( nipij )[i] / sum( psicnipij ) # formula in the text
        # ( 1- zz[,1] ) * ( 1 - zz[,2] ) * rowSums( nipij )[i] / sum( psicnipij ) # wrong index in the text
        ( 1- zz[,1] ) * ( 1 - zz[,2] ) * eta[i] / sum( psicnipij ) # my solution
    }

    # Calculate jacobian for estimating the variance of psi parameters
    # not sure about 5.77
    jpsi <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jpsi[i] <-
        - sum( ww * rowSums( y1y2[,i,] ) / psi[i]^2 ) -
        # sum( ww  * ( 1 - zz[,1] ) * rowSums( sweep( yy[,,2] , 2 , ( nipij[i,] / colSums( sweep( nipij , 1 , (1 - psi)^2 , "*" ) ) ) , "*" ) ) ) - # 2nd attempt
        # sum( ww  * ( 1 - zz[,1] ) * rowSums( sweep( yy[,,2] , 2 , ( nipij[i,] / colSums( psicnipij )^2 ) , "*" ) ) ) - # 3rd attempt
        sum( ww  * ( 1 - zz[,1] ) * rowSums( sweep( yy[,,2] , 2 , ( nipij[i,]^2 / colSums( psicnipij )^2 ) , "*" ) ) ) - # 4th attempt
        sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * eta[i]^2 / sum( psicnipij )^2 ) # my solution
    }

    # divide u_psi by the jacobian
    u_psi <- sweep( u_psi , 2 , jpsi , "/" )

    # calculate variance of psi
    psi_var <- survey::svyrecvar( sweep( u_psi , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    psi_var <- diag( psi_var )

    ### eta

    # Calculate scores for estimating the variance of eta parameters
    u_eta <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_eta[,i] <-
        yy[,i,1] * rowSums( yy[,,2] ) / eta[i] +
        yy[,i,1] * ( 1 - zz[,2] ) / eta[i] +
        rowSums( sweep( yy[,,2] * ( 1 - zz[,1] ) , 2 , psicpij[i,] / colSums( psicnipij ) , "*" ) ) +
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * ( 1 - psi[i] ) / sum( psicni )
    }

    # Calculate jacobian for estimating the variance of eta parameters
    jeta <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jeta[i] <-
        - sum( ww * yy[,i,1] * rowSums( yy[,,2] ) ) / eta[i]^2 -
        sum( ww * yy[,i,1] * ( 1 - zz[,2] ) ) / eta[i]^2 -
        sum( ww * rowSums( sweep( yy[,,2] * ( 1 - zz[,1] ) , 2 , psicpij[i,]^2 / colSums( psicnipij )^2 , "*" ) ) ) -
        sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * ( 1 - psi[i] )^2 / sum( psicni )^2 )
    }

    # divide u_eta by the jacobian
    u_eta <- sweep( u_eta , 2 , jeta , "/" )

    # calculate variance of eta
    eta_var <- survey::svyrecvar( sweep( u_eta , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    eta_var <- diag( eta_var )

    ### pij

    # Calculate scores for estimating the variance of pij parameters
    u_pij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_pij[,i,j] <- ( y1y2[,i,j] / pij[i,j] ) +
        ( yy[,i,1] * ( 1 - zz[,2] ) ) +
        ( yy[,j,2] * ( 1 - zz[,1] ) ) * ( eta[i] / colSums( nipij )[j] ) +
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * eta[i]
    }

    # Calculate jacobian for estimating the variance of pij parameters
    jpij <- array( 0 , dim = c( nrow(bigNij) , ncol(bigNij) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      jpij[i,j] <-
        - sum( ww * y1y2[,i,j] ) / pij[i,j]^2 -
        sum( ww * yy[,i,1] * ( 1 - zz[,2] ) ) / colSums( pij )[j]^2 -
        sum( ww * yy[,j,2] * ( 1 - zz[,1] ) * ( 1 - psi[i] )^2 * eta[i]^2 / sum( psicnipij )^2 )
    }

    # divide u_pij by the jacobian
    u_pij <- sweep( u_pij , 2:3 , jpij , "/" )

    # create matrix of linearized variables
    lin_pij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_pij[,z,] ) )

    # adjust for structural zeros in transition matrix
    if ( !is.null( pij.zero ) ) {
      lin_pij[ , c( t( pij.zero ) ) ] <- 0
    }

    # calculate variance of u_pij sum
    pij_var <- diag( survey::svyrecvar( ww * lin_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    pij_var <- matrix( pij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    # create matrix of linearized variables
    lin_muij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_muij[,z,] ) )

    # adjust for structural zeros in transition matrix
    if ( !is.null( pij.zero ) ) {
      lin_muij[ , c( t( pij.zero ) ) ] <- 0
    }

    # calculate variance of muij
    muij_var <- diag( survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    muij_var <- matrix( muij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### combine results

    # return final estimates
    resvar <-
      list(
        "psi" = psi_var ,
        "rhoRR" = rhoRR_var ,
        "rhoMM" = rhoMM_var ,
        "eta" = eta_var ,
        "pij" = pij_var ,
        "muij" = muij_var )

    # return list of results
    return( resvar )

  } , C = {

    # calculate auxiliary stats
    nipij <- sweep( pij , 1 , eta , "*" )
    rhoni <- rhoMM * eta
    rhocni <- ( 1 - rhoMM ) * eta
    rhonipij <- sweep( pij , 1 , rhoni , "*" )
    rhocnipij <- sweep( pij , 1 , rhocni , "*" )
    rhocpij <- sweep( pij , 1 , 1 - rhoMM , "*" )

    ### psi

    # special variables
    a1 <- ( sum( bigCj ) + bigM ) / N^2
    a2 <- ( sum( bigNij ) + sum( bigRi ) ) / N^2

    # linearized variable e_psi
    e_psi <- a1 * ( 2 - zz[,2] ) + a2 * ( 1 - zz[,1] ) * ( 2 - zz[,2] )

    # calculate variance
    psi_var <- survey::svyrecvar( ww * e_psi , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### rhoRR

    # Calculate scores for estimating the variance of rhoRR parameters
    u_rhoRR <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_rhoRR[,i] <- yy[,i,1] * rowSums( yy[,,2] ) / rhoRR[i] - yy[,i,1] * ( 1 - zz[,2] ) / ( 1 - rhoRR[i] )
    }

    # Calculate jacobian for estimating the variance of rhoRR parameters
    jrhoRR <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jrhoRR[i] <-
        sum( ww * yy[,i,1] * rowSums( yy[,,2] ) / rhoRR[i]^2 ) + sum( ww * yy[,i,1] * ( 1 - zz[,2] ) / ( 1 - rhoRR[i] )^2 )
    }

    # divide u_rhoRR by the jacobian
    u_rhoRR <- sweep( u_rhoRR , 2 , jrhoRR , "/" )

    # calculate variance of rhoRR
    rhoRR_var <- survey::svyrecvar( sweep( u_rhoRR , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    rhoRR_var <- diag( rhoRR_var )

    ### rhoMM

    # Calculate scores for estimating the variance of rhoMM parameters
    u_rhoMM <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_rhoMM[,i] <-
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * eta[i] / sum( rhoni ) -
        rowSums( sweep( yy[,,2] * ( 1 - zz[,1] ) , 2 , nipij[i,] / colSums( rhocnipij ) , "*" ) )
    }

    # Calculate jacobian for estimating the variance of rhoMM parameters
    jrhoMM <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jrhoMM[i] <-
        - sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * eta[i]^2 / sum( rhoni )^2 ) +
        sum( ww * rowSums( sweep( yy[,,2] * ( 1 - zz[,1] ) , 2 , nipij[i,]^2 / colSums( rhocnipij )^2 , "*" ) ) )
    }

    # divide u_rhoMM by the jacobian
    u_rhoMM <- sweep( u_rhoMM , 2 , jrhoMM , "/" )

    # calculate variance of rhoMM
    rhoMM_var <- survey::svyrecvar( sweep( u_rhoMM , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    rhoMM_var <- diag( rhoMM_var )

    ### eta

    # Calculate scores for estimating the variance of eta parameters
    u_eta <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_eta[,i] <-
        rowSums( y1y2[,i,] ) / eta[i] +
        yy[,i,1] * ( 1 - zz[,2] ) / eta[i] +
        rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , rhocpij[i,] / colSums( rhocnipij ) , "*" ) ) +
        ( 1- zz[,1] ) * ( 1 - zz[,2] ) * ( rhoMM[i] / sum( rhoMM ) ) - 1
    }

    # Calculate jacobian for estimating the variance of eta parameters
    jeta <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jeta[i] <-
        - sum( ww * y1y2[,i,] ) / eta[i]^2 - sum( ww * yy[,i,1] * ( 1 - zz[,2] ) ) / eta[i]^2 -
        sum( ww * ( 1 - zz[,1] ) * sweep( yy[,,2] , 2 , sweep( rhocpij^2 , 2 , colSums( rhocnipij )^2 , "/" )[i,] , "*" ) ) -
        sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * rhoni[i]^2 / sum( rhoni )^2 )
    }

    # divide u_eta by the jacobian
    u_eta <- sweep( u_eta , 2 , jeta , "/" )

    # calculate variance of eta
    eta_var <- survey::svyrecvar( sweep( u_eta , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    eta_var <- diag( eta_var )

    ### pij

    # Calculate scores for estimating the variance of pij parameters
    u_pij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_pij[,i,j] <- ( y1y2[,i,j] / pij[i,j] ) +
        ( yy[,i,1] * ( 1 - zz[,2] ) ) +
        ( yy[,j,2] * ( 1 - zz[,1] ) ) * ( rhocni )[i] / colSums( rhocnipij )[j] +
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * ( rhocni )[i] / sum( rhocni )
    }

    # Calculate jacobian for estimating the variance of pij parameters
    jpij <- array( 0 , dim = c( nrow(bigNij) , ncol(bigNij) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      jpij[i,j] <-
        - sum( ww * y1y2[,i,j] ) / pij[i,j]^2 -
        sum( ww * yy[,i,1] * ( 1 - zz[,2] ) ) / colSums( pij )[j]^2 -
        sum( ww * ( 1 - zz[,1] ) * yy[,j,2] ) * ( rhocni )[i]^2 / colSums( rhocnipij )[j]^2
    }

    # divide u_pij by the jacobian
    u_pij <- sweep( u_pij , 2:3 , jpij , "/" )

    # create matrix of linearized variables
    lin_pij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_pij[,z,] ) )

    # adjust for structural zeros in transition matrix
    if ( !is.null( pij.zero ) ) {
      lin_pij[ , c( t( pij.zero ) ) ] <- 0
    }

    # calculate variance of u_pij sum
    pij_var <- diag( survey::svyrecvar( ww * lin_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    pij_var <- matrix( pij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    # create matrix of linearized variables
    lin_muij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_muij[,z,] ) )

    # adjust for structural zeros in transition matrix
    if ( !is.null( pij.zero ) ) {
      lin_muij[ , c( t( pij.zero ) ) ] <- 0
    }

    # calculate variance of muij
    muij_var <- diag( survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    muij_var <- matrix( muij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### combine results

    # return final estimates
    resvar <-
      list(
        "psi" = psi_var ,
        "rhoRR" = rhoRR_var ,
        "rhoMM" = rhoMM_var ,
        "eta" = eta_var ,
        "pij" = pij_var ,
        "muij" = muij_var )

    # return list of results
    return( resvar )

  } ,
  D = {

    # calculate auxiliary stats
    nipij <- sweep( pij , 1 , eta , "*" )
    rhoRRpij <- sweep( pij , 2 , rhoRR , "*" )
    rhoRRcpij <- sweep( pij , 2 , 1 - rhoRR , "*" )
    rhoMMnipij <- sweep( nipij , 2 , rhoMM , "*" )
    rhoMMpij <- sweep( pij , 2 , rhoMM , "*" )
    rhoMMcpij <- sweep( 1 - pij , 2 , rhoMM , "*" )

    ### psi

    # special variables
    a1 <- ( sum( bigCj ) + bigM ) / N^2
    a2 <- ( sum( bigNij ) + sum( bigRi ) ) / N^2

    # linearized variable e_psi
    e_psi <- a1 * ( 2 - zz[,2] ) + a2 * ( 1 - zz[,1] ) * ( 2 - zz[,2] )

    # calculate variance
    psi_var <- survey::svyrecvar( ww * e_psi , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### rhoRR

    # Calculate scores for estimating the variance of rhoRR parameters
    u_rhoRR <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      u_rhoRR[,j] <-
        yy[,j,2] * rowSums( yy[,,1] ) / rhoRR[j] -
        rowSums( sweep( ( 1 - zz[,2] ) * yy[,,1] , 2 ,  pij[,j] / rowSums( rhoMMcpij ) , "*" ) )
    }

    # Calculate jacobian for estimating the variance of rhoRR parameters
    jrhoRR <- vector( "numeric" , length = nrow( bigNij ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      jrhoRR[j] <-
        - sum( ww * yy[,j,2] * rowSums( yy[,,1] ) ) / rhoRR[j]^2 +
        sum( ww * sweep( yy[,,1] * ( 1 - zz[,2] ) , 2 , pij[,j]^2 / rowSums( rhoMMcpij )^2 , "*" ) )
    }

    # divide u_rhoRR by the jacobian
    u_rhoRR <- sweep( u_rhoRR , 2 , jrhoRR , "/" )

    # calculate variance of rhoRR
    rhoRR_var <- survey::svyrecvar( sweep( u_rhoRR , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    rhoRR_var <- diag( rhoRR_var )

    ### rhoMM

    # Calculate scores for estimating the variance of rhoMM parameters
    u_rhoMM <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      u_rhoMM[,j] <-
        - yy[,j,2] * ( 1 - zz[,1] ) / ( 1 - rhoMM[j] ) +
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * sum( nipij[,j] ) / sum( rhoMMnipij )
    }

    # Calculate jacobian for estimating the variance of rhoMM parameters
    jrhoMM <- vector( "numeric" , length = nrow( bigNij ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      jrhoMM[j] <-
        - sum( ww * yy[,j,2] * ( 1 - zz[,1] ) ) / ( 1 - rhoMM[j] )^2 -
        sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * sum( nipij[,j] )^2 / sum( rhoMMnipij )^2 )
    }

    # divide u_rhoMM by the jacobian
    u_rhoMM <- sweep( u_rhoMM , 2 , jrhoMM , "/" )

    # calculate variance of rhoMM
    rhoMM_var <- survey::svyrecvar( sweep( u_rhoMM , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    rhoMM_var <- diag( rhoMM_var )

    ### eta

    # Calculate scores for estimating the variance of eta parameters
    u_eta <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_eta[,i] <-
        rowSums( y1y2[,i,] ) / eta[i] +
        yy[,i,1] * ( 1 - zz[,2] ) / eta[i] +
        rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , pij[i,] / colSums( nipij ) , "*" ) ) +
        ( 1- zz[,1] ) * ( 1 - zz[,2] ) * ( rowSums( rhoMMpij )[i] / sum( rhoMMnipij ) ) - 1
    }

    # Calculate jacobian for estimating the variance of eta parameters
    jeta <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jeta[i] <-
        - sum( ww * y1y2[,i,] ) / eta[i]^2 - sum( ww * yy[,i,1] * ( 1 - zz[,2] ) ) / eta[i]^2 -
        sum( ww * ( 1 - zz[,1] ) * sweep( yy[,,2] , 2 , sweep( pij^2 , 2 , colSums( nipij )^2 , "/" )[i,] , "*" ) ) -
        sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * rhoMMpij[i]^2 / sum( rhoMMnipij )^2 )
    }

    # divide u_eta by the jacobian
    u_eta <- sweep( u_eta , 2 , jeta , "/" )

    # calculate variance of eta
    eta_var <- survey::svyrecvar( sweep( u_eta , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    eta_var <- diag( eta_var )

    ### pij

    # Calculate scores for estimating the variance of pij parameters
    # lambda2 <- -( rowSums( bigNij ) + bigRi + rowSums( sweep( nipij , 2 , bigCj / colSums( nipij ) , "*" ) ) + bigM * rowSums( rhoMMnipij ) / sum( rhoMMnipij ) ) / N
    u_pij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_pij[,i,j] <- ( y1y2[,i,j] / pij[i,j] ) +
        yy[,i,1] * ( 1 - zz[,2] ) * ( 1 - rhoRR[j] ) / rowSums( rhoRRcpij )[i] +
        yy[,j,2] * ( 1 - zz[,1] ) * ( eta )[i] / colSums( nipij )[j] +
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * ( rhoMM[j] * eta[i] ) / sum( rhoMMnipij )
    }

    # Calculate jacobian for estimating the variance of pij parameters
    jpij <- array( 0 , dim = c( nrow(bigNij) , ncol(bigNij) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      jpij[i,j] <-
        - sum( ww * y1y2[,i,j] ) / pij[i,j]^2 -
        sum( ww * yy[,i,1] * ( 1 - zz[,2] ) * ( 1 - rhoRR[j] )^2 ) / rowSums( rhoRRcpij )[i]^2 -
        sum ( ww * yy[,j,2] * ( 1 - zz[,1] ) * ( eta )[i]^2 ) / colSums( nipij )[j]^2 -
        sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * ( rhoMM[j] * eta[i] )^2 / sum( rhoMMnipij )^2 )
    }

    # divide u_pij by the jacobian
    u_pij <- sweep( u_pij , 2:3 , jpij , "/" )

    # create matrix of linearized variables
    lin_pij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_pij[,z,] ) )

    # adjust for structural zeros in transition matrix
    if ( !is.null( pij.zero ) ) {
      lin_pij[ , c( t( pij.zero ) ) ] <- 0
    }

    # calculate variance of u_pij sum
    pij_var <- diag( survey::svyrecvar( ww * lin_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    pij_var <- matrix( pij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    # create matrix of linearized variables
    lin_muij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_muij[,z,] ) )

    # adjust for structural zeros in transition matrix
    if ( !is.null( pij.zero ) ) {
      lin_muij[ , c( t( pij.zero ) ) ] <- 0
    }

    # calculate variance of muij
    muij_var <- diag( survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    muij_var <- matrix( muij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### combine results

    # return final estimates
    resvar <-
      list(
        "psi" = psi_var ,
        "rhoRR" = rhoRR_var ,
        "rhoMM" = rhoMM_var ,
        "eta" = eta_var ,
        "pij" = pij_var ,
        "muij" = muij_var )

    # return list of results
    return( resvar )

  } )

  # return variance
  mvar

}
