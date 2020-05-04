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
    muij_var <- survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    des_vcov <- muij_var # full variance-covariance matrix
    muij_var <- matrix( diag( muij_var ) , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    # calculate variance under SRS with replacement
    # for Kish DEff
    nobs <- sum( ww > 0 )
    srs_vcov <- svyvar( ( ww > 0 ) * lin_muij , design, na.rm = TRUE ) * sum( ww )^2 / nobs
    srs_vcov <- matrix( diag( srs_vcov ) , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### combine results

    # return final estimates
    resvar <-
      list(
        "psi" = psi_var ,
        "rhoRR" = rhoRR_var ,
        "rhoMM" = rhoMM_var ,
        "eta" = eta_var ,
        "pij" = pij_var ,
        "muij" = muij_var ,
        "vcov_full" = des_vcov ,
        "vcov_srs"  = srs_vcov )

  } , A = {

    # calculate auxiliary stats
    nipij <- sweep( pij , 1 , eta , "*" )

    ### psi

    # Calculate scores for estimating the variance of psi parameters
    u_psi <-
      rowSums( apply( y1y2 , 1:2 , sum ) ) / psi +
      rowSums( yy[,,1] * (1 - zz[,2] ) ) / psi -
      rowSums( yy[,,2] * (1 - zz[,1] ) ) / ( 1 - psi ) -
      ( 1- zz[,1] ) * ( 1 - zz[,2] ) / ( 1 - psi )

    # Calculate jacobian for estimating the variance of psi parameters
    jpsi <- - sum( bigNij ) / psi^2 - sum( bigRi ) / psi^2 - sum( bigCj ) / ( 1 - psi )^2 - bigM / ( 1 - psi )^2

    # divide u_psi by the jacobian
    u_psi <- u_psi / jpsi

    # calculate variance of psi
    psi_var <- survey::svyrecvar( ww * u_psi , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### rhoRR

    # Calculate scores for estimating the variance of rhoRR parameters
    u_rhoRR <- rowSums( apply( y1y2 , 1:2 , sum ) ) / rhoRR - rowSums( yy[,,1] * ( 1 - zz[,2]) ) / ( 1 - rhoRR )

    # Calculate jacobian for estimating the variance of rhoRR parameters
    jrhoRR <- - sum( bigNij ) / rhoRR^2 - sum( bigRi ) / ( 1 - rhoRR )^2

    # divide u_rhoRR by the jacobian
    u_rhoRR <- u_rhoRR / jrhoRR

    # calculate variance of rhoRR
    rhoRR_var <- survey::svyrecvar( ww * u_rhoRR , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### rhoMM

    # Calculate scores for estimating the variance of rhoMM parameters
    u_rhoMM <- ( 1 - zz[,1] ) * ( 1 - zz[,2] ) / rhoMM - rowSums( yy[,,2] * ( 1 - zz[,1] ) ) / ( 1 - rhoMM )

    # Calculate jacobian for estimating the variance of rhoMM parameters
    jrhoMM <- - bigM / rhoMM^2 - sum( bigCj ) / ( 1 - rhoMM )^2

    # divide u_rhoMM by the jacobian
    u_rhoMM <- u_rhoMM / jrhoMM

    # calculate variance of rhoMM
    rhoMM_var <- survey::svyrecvar( ww * u_rhoMM , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### eta

    # Calculate scores for estimating the variance of eta parameters
    u_eta <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_eta[,i] <-
        rowSums( y1y2[,i,] ) / eta[i] +
        yy[,i,1] * ( 1 - zz[,2] ) / eta[i] +
        rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , pij[i,] / colSums( nipij ) , "*" ) ) +
        ( 1- zz[,1] ) * ( 1 - zz[,2] )
    }

    # Calculate jacobian for estimating the variance of eta parameters
    jeta <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jeta[i] <- - sum( bigNij[i,] ) / eta[i]^2 - sum( bigCj * pij[i,]^2 / colSums( nipij )[i]^2 )
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
      jpij[i,j] <- - bigNij[i,j] / pij[i,j]^2 - bigCj[j] * eta[i]^2 / colSums( nipij )[j]^2
      # jpij[i,j] <- - bigNij[i,j] / pij[i,j]^2 - bigRi[i] - bigCj[j] * eta[i]^2 / colSums( nipij )[j]^2
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
    muij_var <- survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    des_vcov <- muij_var # full variance-covariance matrix
    muij_var <- matrix( diag( muij_var ) , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    # calculate variance under SRS with replacement
    # for Kish DEff
    nobs <- sum( ww > 0 )
    srs_vcov <- svyvar( ( ww > 0 ) * lin_muij , design, na.rm = TRUE ) * sum( ww )^2 / nobs
    srs_vcov <- matrix( diag( srs_vcov ) , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### combine results

    # return final estimates
    resvar <-
      list(
        "psi" = psi_var ,
        "rhoRR" = rhoRR_var ,
        "rhoMM" = rhoMM_var ,
        "eta" = eta_var ,
        "pij" = pij_var ,
        "muij" = muij_var ,
        "vcov_full" = des_vcov ,
        "vcov_srs"  = srs_vcov )

  } , B ={

    # calculate auxiliary stats
    nipij <- sweep( pij , 1 , eta , "*" )
    psicni <- ( 1 - psi ) * eta
    psicnipij <- sweep( pij , 1 , psicni , "*" )
    psicpij <- sweep( pij , 1 , 1 - psi , "*" )

    ### psi

    # Calculate scores for estimating the variance of psi parameters
    u_psi <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_psi[,i] <-
        yy[,i,1] * rowSums( yy[,,2] ) / psi[i] +
        yy[,i,1] * ( 1 - zz[,2] ) / psi[i] -
        rowSums( sweep( ( 1 - zz[,1] ) * yy[,,2] , 2 , nipij[i,] / colSums( psicnipij ) , "*" ) ) -
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * eta[i] / sum( psicni )
    }

    # Calculate jacobian for estimating the variance of psi parameters
    jpsi <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jpsi[i] <-
        - sum( bigNij[i,] ) / psi[i]^2 -
        bigRi[i] / psi[i]^2 -
        sum( bigCj * nipij[i,]^2 / colSums( psicnipij )^2 ) -
        bigM * eta[i]^2 / sum( psicni )^2
    }

    # divide u_psi by the jacobian
    u_psi <- sweep( u_psi , 2 , jpsi , "/" )

    # calculate variance of psi
    psi_var <- survey::svyrecvar( sweep( u_psi , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    psi_var <- diag( psi_var )

    ### rhoRR

    # Calculate scores for estimating the variance of rhoRR parameters
    u_rhoRR <- rowSums( apply( y1y2 , 1:2 , sum ) ) / rhoRR - rowSums( yy[,,1] * ( 1 - zz[,2]) ) / ( 1 - rhoRR )

    # Calculate jacobian for estimating the variance of rhoRR parameters
    jrhoRR <- - sum( bigNij ) / rhoRR^2 - sum( bigRi ) / ( 1 - rhoRR )^2

    # divide u_rhoRR by the jacobian
    u_rhoRR <- u_rhoRR / jrhoRR

    # calculate variance of rhoRR
    rhoRR_var <- survey::svyrecvar( ww * u_rhoRR , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### rhoMM

    # Calculate scores for estimating the variance of rhoMM parameters
    u_rhoMM <- ( 1 - zz[,1] ) * ( 1 - zz[,2] ) / rhoMM - rowSums( yy[,,2] * ( 1 - zz[,1] ) ) / ( 1 - rhoMM )

    # Calculate jacobian for estimating the variance of rhoMM parameters
    jrhoMM <- - bigM / rhoMM^2 - sum( bigCj ) / ( 1 - rhoMM )^2

    # divide u_rhoMM by the jacobian
    u_rhoMM <- u_rhoMM / jrhoMM

    # calculate variance of rhoMM
    rhoMM_var <- survey::svyrecvar( ww * u_rhoMM , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### eta

    # Calculate scores for estimating the variance of eta parameters
    u_eta <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_eta[,i] <-
        yy[,i,1] * rowSums( yy[,,2] ) / eta[i] +
        yy[,i,1] * ( 1 - zz[,2] ) / eta[i] +
        rowSums( sweep( ( 1 - zz[,1] ) * yy[,,2] , 2 , psicpij[i,] / colSums( psicnipij ) , "*" ) ) +
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * ( 1 - psi[i] ) / sum( psicni )
    }

    # Calculate jacobian for estimating the variance of eta parameters
    jeta <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jeta[i] <-
        - sum( bigNij[i,] ) / eta[i]^2 -
        bigRi[i] / eta[i]^2 -
        sum( bigCj * psicpij[i,]^2 / colSums( psicnipij )^2 ) -
        bigM * ( 1- psi[i] )^2 / sum( psicni )^2
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
        ( yy[,j,2] * ( 1 - zz[,1] ) ) * ( psicni[i] / colSums( psicnipij )[j] ) +
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * ( psicni[i] / sum( psicnipij ) )
    }

    # Calculate jacobian for estimating the variance of pij parameters
    jpij <- array( 0 , dim = c( nrow(bigNij) , ncol(bigNij) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      jpij[i,j] <-
        - bigNij[i,j] / pij[i,j]^2 -
        bigRi[i] -
        bigCj[j] * psicni[i]^2 / colSums( psicnipij )[j]^2 -
        bigM * ( psicni[i]^2 ) / sum( psicni )^2
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
    muij_var <- survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    des_vcov <- muij_var # full variance-covariance matrix
    muij_var <- matrix( diag( muij_var ) , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    # calculate variance under SRS with replacement
    # for Kish DEff
    nobs <- sum( ww > 0 )
    srs_vcov <- svyvar( ( ww > 0 ) * lin_muij , design, na.rm = TRUE ) * sum( ww )^2 / nobs
    srs_vcov <- matrix( diag( srs_vcov ) , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### combine results

    # return final estimates
    resvar <-
      list(
        "psi" = psi_var ,
        "rhoRR" = rhoRR_var ,
        "rhoMM" = rhoMM_var ,
        "eta" = eta_var ,
        "pij" = pij_var ,
        "muij" = muij_var ,
        "vcov_full" = des_vcov ,
        "vcov_srs"  = srs_vcov )

  } , C = {

    # calculate auxiliary stats
    nipij <- sweep( pij , 1 , eta , "*" )
    rhoni <- rhoMM * eta
    rhocni <- ( 1 - rhoMM ) * eta
    rhonipij <- sweep( pij , 1 , rhoni , "*" )
    rhocnipij <- sweep( pij , 1 , rhocni , "*" )
    rhocpij <- sweep( pij , 1 , 1 - rhoMM , "*" )

    ### psi

    # Calculate scores for estimating the variance of psi parameters
    u_psi <-
      rowSums( apply( y1y2 , 1:2 , sum ) ) / psi +
      rowSums( yy[,,1] * (1 - zz[,2] ) ) / psi -
      rowSums( yy[,,2] * (1 - zz[,1] ) ) / ( 1 - psi ) -
      ( 1- zz[,1] ) * ( 1 - zz[,2] ) / ( 1 - psi )

    # Calculate jacobian for estimating the variance of psi parameters
    jpsi <- - sum( bigNij ) / psi^2 - sum( bigRi ) / psi^2 - sum( bigCj ) / ( 1 - psi )^2 - bigM / ( 1 - psi )^2

    # divide u_psi by the jacobian
    u_psi <- u_psi / jpsi

    # calculate variance of psi
    psi_var <- survey::svyrecvar( ww * u_psi , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### rhoRR

    # Calculate scores for estimating the variance of rhoRR parameters
    u_rhoRR <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_rhoRR[,i] <- yy[,i,1] * rowSums( yy[,,2] ) / rhoRR[i] - yy[,i,1] * ( 1 - zz[,2] ) / ( 1 - rhoRR[i] )
    }

    # Calculate jacobian for estimating the variance of rhoRR parameters
    jrhoRR <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jrhoRR[i] <- - sum( bigNij[i,] ) / rhoRR[i]^2 - bigRi[i] / ( 1 - rhoRR[i] )^2
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
         bigM * eta[i]^2 / sum( rhoni )^2 -
        sum( bigCj * nipij[i,]^2 / colSums( rhocnipij )^2 )
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
        yy[,i,1] * rowSums( yy[,,2] ) / eta[i] +
        yy[,i,1] * ( 1 - zz[,2] ) / eta[i] +
        rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , rhocpij[i,] / colSums( rhocnipij ) , "*" ) ) +
        ( 1- zz[,1] ) * ( 1 - zz[,2] ) * ( rhoMM[i] / sum( rhoni ) ) - 1
    }

    # Calculate jacobian for estimating the variance of eta parameters
    jeta <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jeta[i] <-
        - sum( bigNij[i,] ) / eta[i]^2 -
        bigRi[i] / eta[i]^2 -
        sum( bigCj * rhocpij[i,]^2 / colSums( rhocnipij )^2 ) -
        bigM * rhoMM[i]^2 / sum( rhoni )^2
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
        - bigNij[i,j] / pij[i,j]^2 -
        bigRi[i] -
        bigCj[j] * ( rhocni )[i]^2 / colSums( rhocnipij )[j]^2 -
        bigM * rhoni[i]^2 / sum( rhoni )^2
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
    muij_var <- survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    des_vcov <- muij_var # full variance-covariance matrix
    muij_var <- matrix( diag( muij_var ) , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    # calculate variance under SRS with replacement
    # for Kish DEff
    nobs <- sum( ww > 0 )
    srs_vcov <- svyvar( ( ww > 0 ) * lin_muij , design, na.rm = TRUE ) * sum( ww )^2 / nobs
    srs_vcov <- matrix( diag( srs_vcov ) , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### combine results

    # return final estimates
    resvar <-
      list(
        "psi" = psi_var ,
        "rhoRR" = rhoRR_var ,
        "rhoMM" = rhoMM_var ,
        "eta" = eta_var ,
        "pij" = pij_var ,
        "muij" = muij_var ,
        "vcov_full" = des_vcov ,
        "vcov_srs"  = srs_vcov )

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

    # Calculate scores for estimating the variance of psi parameters
    u_psi <-
      rowSums( apply( y1y2 , 1:2 , sum ) ) / psi +
      rowSums( yy[,,1] * (1 - zz[,2] ) ) / psi -
      rowSums( yy[,,2] * (1 - zz[,1] ) ) / ( 1 - psi ) -
      ( 1- zz[,1] ) * ( 1 - zz[,2] ) / ( 1 - psi )

    # Calculate jacobian for estimating the variance of psi parameters
    jpsi <- - sum( bigNij ) / psi^2 - sum( bigRi ) / psi^2 - sum( bigCj ) / ( 1 - psi )^2 - bigM / ( 1 - psi )^2

    # divide u_psi by the jacobian
    u_psi <- u_psi / jpsi

    # calculate variance of psi
    psi_var <- survey::svyrecvar( ww * u_psi , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### rhoRR

    # Calculate scores for estimating the variance of rhoRR parameters
    u_rhoRR <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      u_rhoRR[,j] <-
        yy[,j,2] * rowSums( yy[,,1] ) / rhoRR[j] -
        rowSums( sweep( ( 1 - zz[,2] ) * yy[,,1] , 2 ,  pij[,j] / rowSums( rhoRRcpij ) , "*" ) )
    }

    # Calculate jacobian for estimating the variance of rhoRR parameters
    jrhoRR <- vector( "numeric" , length = nrow( bigNij ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      jrhoRR[j] <-
        - sum( bigNij[,j] ) / rhoRR[j]^2 -
        sum( bigRi * pij[,j]^2 / rowSums( rhoRRcpij )^2 )
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
         - bigCj[j] / ( 1 - rhoMM[j] )^2 - bigM * sum( nipij[,j] )^2 / sum( rhoMMnipij )^2
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
        yy[,i,1] * rowSums( yy[,,2] ) / eta[i] +
        yy[,i,1] * ( 1 - zz[,2] ) / eta[i] +
        rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , pij[i,] / colSums( nipij ) , "*" ) ) +
        ( 1- zz[,1] ) * ( 1 - zz[,2] ) * ( rowSums( rhoMMpij )[i] / sum( rhoMMnipij ) ) - 1
    }

    # Calculate jacobian for estimating the variance of eta parameters
    jeta <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jeta[i] <-
        - sum( bigNij[i,] ) / eta[i]^2 -
        bigRi[i] / eta[i]^2 -
        sum( bigCj * pij[i,]^2 / colSums( nipij )^2 ) -
        bigM * rowSums( rhoMMpij )[i]^2 / sum( rhoMMnipij )^2
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
        - bigNij[i,j] / pij[i,j]^2 -
        bigRi[i] * ( 1 - rhoRR[j] )^2 / rowSums( rhoRRcpij )[i]^2 -
        bigCj[j] * eta[i]^2 / colSums( nipij )[j]^2
        bigM * ( rhoMM[j] * eta[i] )^2 / sum( rhoMMnipij )^2
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
    muij_var <- survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    des_vcov <- muij_var # full variance-covariance matrix
    muij_var <- matrix( diag( muij_var ) , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    # calculate variance under SRS with replacement
    # for Kish DEff
    nobs <- sum( ww > 0 )
    srs_vcov <- svyvar( ( ww > 0 ) * lin_muij , design, na.rm = TRUE ) * sum( ww )^2 / nobs
    srs_vcov <- matrix( diag( srs_vcov ) , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### combine results

    # return final estimates
    resvar <-
      list(
        "psi" = psi_var ,
        "rhoRR" = rhoRR_var ,
        "rhoMM" = rhoMM_var ,
        "eta" = eta_var ,
        "pij" = pij_var ,
        "muij" = muij_var ,
        "vcov_full" = des_vcov ,
        "vcov_srs"  = srs_vcov )

  } )

  # calculate 1st rao-scott correction
  if ( res$model != "MCAR" ) {

    # # F-Distribution        # not working
    # # DeltaMat <- Matrix::solve( mvar$vcov_full ) %*% mvar$vcov_srs
    # d0 <- sum( diag( DeltaMat ) )^2 / ( sum( diag( DeltaMat %*% DeltaMat ) ) )
    # nu <- length( unique( design$cluster[, 1] ) ) - length( unique( design$strata[,1] ) )
    # pearson <- chisq.test( res$estimated.counts , correct = FALSE )
    # pearson$statistic <- pearson$statistic / sum( diag( DeltaMat ) )
    # pearson$p.value <- pf( pearson$statistic, d0, d0 * nu, lower.tail = FALSE )
    # attr( pearson$statistic, "names" ) <- "F"
    # pearson$parameter <- c( ndf = d0 , ddf = d0 * nu)
    # pearson$method <- "Pearson's X^2: 2nd Rao & Scott adjustment"
    # pearson$data.name <- "estimated counts"

    # Chi-2 Distribution
    DeltaMat <- Matrix::solve( mvar$muij ) %*% mvar$vcov_srs
    pearson <- chisq.test( res$estimated.counts , correct = FALSE )
    pearson$statistic <- pearson$statistic / sum( diag( DeltaMat ) )
    pearson$p.value <- pchisq(pearson$statistic/mean(diag(DeltaMat)), df = NCOL(DeltaMat) , lower.tail = FALSE)
    pearson$parameter <- c(df = NCOL(DeltaMat))
    pearson$method <- "Pearson's X^2: 1st Rao & Scott adjustment"
    pearson$data.name <- "estimated counts"

  }

  # add results
  mvar[["adj.chisq"]] <- if ( res$model == "MCAR" ) NULL else pearson

  # return variance
  mvar

}
