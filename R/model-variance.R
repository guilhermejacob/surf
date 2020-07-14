# function for model variance calculation
ipf_variance <- function( xx , ww , res , design , rp.variance = TRUE ) {

  # # load objects
  # for ( this_obj in names( res ) ) assign( this_obj , res[[this_obj]] )

  # load objects
  bigNij <- res[["bigNij"]]
  bigRi <- res[["bigRi"]]
  bigCj <- res[["bigCj"]]
  bigM <- res[["bigM"]]
  psi <- res[["psi"]]
  rho <- res[["rho"]]
  tau <- res[["tau"]]
  eta <- res[["eta"]]
  pij <- res[["pij"]]
  muij <- res[["muij"]]
  pij.zero <- res[["pij.zero"]]
  N <- res[["N"]]
  pearson <- res[["unadj.chisq"]]
  ll <- res[["ll"]]

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
  mlin <- switch( res$model , MCAR = {

    ### psi

    psi_var <- NA

    ### rho

    rho_var <- NA

    ### tau

    tau_var <- NA

    # aux stats
    nipij <- sweep( pij , 1 , eta , "*" )

    ### eta

    # Calculate scores for estimating the variance of eta parameters
    u_eta <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_eta[,i] <-
        ( sum( bigNij ) * rowSums( y1y2[,i,] ) - rowSums( bigNij )[i] * rowSums( y1y2[,i,] ) * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) ) / sum( bigNij )^2
    }

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

    ### gamma

    # calculate linearized variables
    u_gamma <- sapply( seq_len( ncol( pij ) ) , function(j) {
      rowSums( sweep( u_pij[,,j] , 2 , eta , "*" ) + sweep( u_eta , 2 , pij[,j] , "*" ) )
    } )

    ### combine results

    # return linearized variables
    reslin <-
      list(
        "psi" = NA ,
        "rho" = NA ,
        "tau" = NA ,
        "eta" = u_eta ,
        "pij" = u_pij ,
        "muij" = u_muij ,
        "gamma" = u_gamma )

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

    ### rho

    # Calculate scores for estimating the variance of rho parameters
    u_rho <- rowSums( apply( y1y2 , 1:2 , sum ) ) / rho - rowSums( yy[,,1] * ( 1 - zz[,2]) ) / ( 1 - rho )

    # Calculate jacobian for estimating the variance of rho parameters
    jrho <- - sum( bigNij ) / rho^2 - sum( bigRi ) / ( 1 - rho )^2

    # divide u_rho by the jacobian
    u_rho <- u_rho / jrho

    ### tau

    # Calculate scores for estimating the variance of tau parameters
    u_tau <- ( 1 - zz[,1] ) * ( 1 - zz[,2] ) / tau - rowSums( yy[,,2] * ( 1 - zz[,1] ) ) / ( 1 - tau )

    # Calculate jacobian for estimating the variance of tau parameters
    jtau <- - bigM / tau^2 - sum( bigCj ) / ( 1 - tau )^2

    # divide u_tau by the jacobian
    u_tau <- u_tau / jtau

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
      jeta[i] <- - sum( bigNij[i,] ) / eta[i]^2 - sum( bigCj * pij[i,]^2 / colSums( nipij )[i]^2 ) - bigM
    }

    # divide u_eta by the jacobian
    u_eta <- sweep( u_eta , 2 , jeta , "/" )

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
      # jpij[i,j] <- - bigNij[i,j] / pij[i,j]^2 - bigCj[j] * eta[i]^2 / colSums( nipij )[j]^2
      jpij[i,j] <- - bigNij[i,j] / pij[i,j]^2 - bigRi[i] - bigCj[j] * eta[i]^2 / colSums( nipij )[j]^2 - bigM * eta[i]^2
    }

    # divide u_pij by the jacobian
    u_pij <- sweep( u_pij , 2:3 , jpij , "/" )

    ### muij

    # linearized variable
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    ### gamma

    # calculate linearized variables
    u_gamma <- sapply( seq_len( ncol( pij ) ) , function(j) {
      rowSums( sweep( u_pij[,,j] , 2 , eta , "*" ) + sweep( u_eta , 2 , pij[,j] , "*" ) )
    } )

    ### combine results

    # return linearized variables
    reslin <-
      list(
        "psi" = u_psi ,
        "rho" = u_rho ,
        "tau" = u_tau ,
        "eta" = u_eta ,
        "pij" = u_pij ,
        "muij" = u_muij ,
        "gamma" = u_gamma )

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

    ### rho

    # Calculate scores for estimating the variance of rho parameters
    u_rho <- rowSums( apply( y1y2 , 1:2 , sum ) ) / rho - rowSums( yy[,,1] * ( 1 - zz[,2]) ) / ( 1 - rho )

    # Calculate jacobian for estimating the variance of rho parameters
    jrho <- - sum( bigNij ) / rho^2 - sum( bigRi ) / ( 1 - rho )^2

    # divide u_rho by the jacobian
    u_rho <- u_rho / jrho

    ### tau

    # Calculate scores for estimating the variance of tau parameters
    u_tau <- ( 1 - zz[,1] ) * ( 1 - zz[,2] ) / tau - rowSums( yy[,,2] * ( 1 - zz[,1] ) ) / ( 1 - tau )

    # Calculate jacobian for estimating the variance of tau parameters
    jtau <- - bigM / tau^2 - sum( bigCj ) / ( 1 - tau )^2

    # divide u_tau by the jacobian
    u_tau <- u_tau / jtau

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

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    ### gamma

    # calculate linearized variables
    u_gamma <- sapply( seq_len( ncol( pij ) ) , function(j) {
      rowSums( sweep( u_pij[,,j] , 2 , eta , "*" ) + sweep( u_eta , 2 , pij[,j] , "*" ) )
    } )

    ### combine results

    # return linearized variables
    reslin <-
      list(
        "psi" = u_psi ,
        "rho" = u_rho ,
        "tau" = u_tau ,
        "eta" = u_eta ,
        "pij" = u_pij ,
        "muij" = u_muij ,
        "gamma" = u_gamma )

  } , C = {

    # calculate auxiliary stats
    nipij <- sweep( pij , 1 , eta , "*" )
    rhoni <- tau * eta
    rhocni <- ( 1 - tau ) * eta
    rhonipij <- sweep( pij , 1 , rhoni , "*" )
    rhocnipij <- sweep( pij , 1 , rhocni , "*" )
    rhocpij <- sweep( pij , 1 , 1 - tau , "*" )

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

    ### rho

    # Calculate scores for estimating the variance of rho parameters
    u_rho <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_rho[,i] <- yy[,i,1] * rowSums( yy[,,2] ) / rho[i] - yy[,i,1] * ( 1 - zz[,2] ) / ( 1 - rho[i] )
    }

    # Calculate jacobian for estimating the variance of rho parameters
    jrho <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jrho[i] <- - sum( bigNij[i,] ) / rho[i]^2 - bigRi[i] / ( 1 - rho[i] )^2
    }

    # divide u_rho by the jacobian
    u_rho <- sweep( u_rho , 2 , jrho , "/" )

    ### tau

    # Calculate scores for estimating the variance of tau parameters
    u_tau <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_tau[,i] <-
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * eta[i] / sum( rhoni ) -
        rowSums( sweep( yy[,,2] * ( 1 - zz[,1] ) , 2 , nipij[i,] / colSums( rhocnipij ) , "*" ) )
    }

    # Calculate jacobian for estimating the variance of tau parameters
    jtau <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jtau[i] <-
        bigM * eta[i]^2 / sum( rhoni )^2 -
        sum( bigCj * nipij[i,]^2 / colSums( rhocnipij )^2 )
    }

    # divide u_tau by the jacobian
    u_tau <- sweep( u_tau , 2 , jtau , "/" )

    ### eta

    # Calculate scores for estimating the variance of eta parameters
    u_eta <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_eta[,i] <-
        yy[,i,1] * rowSums( yy[,,2] ) / eta[i] +
        yy[,i,1] * ( 1 - zz[,2] ) / eta[i] +
        rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , rhocpij[i,] / colSums( rhocnipij ) , "*" ) ) +
        ( 1- zz[,1] ) * ( 1 - zz[,2] ) * ( tau[i] / sum( rhoni ) ) - 1
    }

    # Calculate jacobian for estimating the variance of eta parameters
    jeta <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jeta[i] <-
        - sum( bigNij[i,] ) / eta[i]^2 -
        bigRi[i] / eta[i]^2 -
        sum( bigCj * rhocpij[i,]^2 / colSums( rhocnipij )^2 ) -
        bigM * tau[i]^2 / sum( rhoni )^2
    }

    # divide u_eta by the jacobian
    u_eta <- sweep( u_eta , 2 , jeta , "/" )

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

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    # gamma

    # calculate linearized variables
    u_gamma <- sapply( seq_len( ncol( pij ) ) , function(j) {
      rowSums( sweep( u_pij[,,j] , 2 , eta , "*" ) + sweep( u_eta , 2 , pij[,j] , "*" ) )
    } )

    ### combine results

    # return linearized variables
    reslin <-
      list(
        "psi" = u_psi ,
        "rho" = u_rho ,
        "tau" = u_tau ,
        "eta" = u_eta ,
        "pij" = u_pij ,
        "muij" = u_muij ,
        "gamma" = u_gamma )

  } ,
  D = {

    # calculate auxiliary stats
    nipij <- sweep( pij , 1 , eta , "*" )
    rhopij <- sweep( pij , 2 , rho , "*" )
    rhocpij <- sweep( pij , 2 , 1 - rho , "*" )
    taunipij <- sweep( nipij , 2 , tau , "*" )
    taupij <- sweep( pij , 2 , tau , "*" )
    taucpij <- sweep( 1 - pij , 2 , tau , "*" )

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

    ### rho

    # Calculate scores for estimating the variance of rho parameters
    u_rho <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      u_rho[,j] <-
        yy[,j,2] * rowSums( yy[,,1] ) / rho[j] -
        rowSums( sweep( ( 1 - zz[,2] ) * yy[,,1] , 2 ,  pij[,j] / rowSums( rhocpij ) , "*" ) )
    }

    # Calculate jacobian for estimating the variance of rho parameters
    jrho <- vector( "numeric" , length = nrow( bigNij ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      jrho[j] <-
        - sum( bigNij[,j] ) / rho[j]^2 -
        sum( bigRi * pij[,j]^2 / rowSums( rhocpij )^2 )
    }

    # divide u_rho by the jacobian
    u_rho <- sweep( u_rho , 2 , jrho , "/" )

    ### tau

    # Calculate scores for estimating the variance of tau parameters
    u_tau <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      u_tau[,j] <-
        - yy[,j,2] * ( 1 - zz[,1] ) / ( 1 - tau[j] ) +
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * sum( nipij[,j] ) / sum( taunipij )
    }

    # Calculate jacobian for estimating the variance of tau parameters
    jtau <- vector( "numeric" , length = nrow( bigNij ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      jtau[j] <-
        - bigCj[j] / ( 1 - tau[j] )^2 - bigM * sum( nipij[,j] )^2 / sum( taunipij )^2
    }

    # divide u_tau by the jacobian
    u_tau <- sweep( u_tau , 2 , jtau , "/" )

    ### eta

    # Calculate scores for estimating the variance of eta parameters
    u_eta <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_eta[,i] <-
        yy[,i,1] * rowSums( yy[,,2] ) / eta[i] +
        yy[,i,1] * ( 1 - zz[,2] ) / eta[i] +
        rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , pij[i,] / colSums( nipij ) , "*" ) ) +
        ( 1- zz[,1] ) * ( 1 - zz[,2] ) * ( rowSums( taupij )[i] / sum( taunipij ) ) - 1
    }

    # Calculate jacobian for estimating the variance of eta parameters
    jeta <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jeta[i] <-
        - sum( bigNij[i,] ) / eta[i]^2 -
        bigRi[i] / eta[i]^2 -
        sum( bigCj * pij[i,]^2 / colSums( nipij )^2 ) -
        bigM * rowSums( taupij )[i]^2 / sum( taunipij )^2
    }

    # divide u_eta by the jacobian
    u_eta <- sweep( u_eta , 2 , jeta , "/" )

    ### pij

    # Calculate scores for estimating the variance of pij parameters
    # lambda2 <- -( rowSums( bigNij ) + bigRi + rowSums( sweep( nipij , 2 , bigCj / colSums( nipij ) , "*" ) ) + bigM * rowSums( taunipij ) / sum( taunipij ) ) / N
    u_pij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_pij[,i,j] <- ( y1y2[,i,j] / pij[i,j] ) +
        yy[,i,1] * ( 1 - zz[,2] ) * ( 1 - rho[j] ) / rowSums( rhocpij )[i] +
        yy[,j,2] * ( 1 - zz[,1] ) * ( eta )[i] / colSums( nipij )[j] +
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * ( tau[j] * eta[i] ) / sum( taunipij )
    }

    # Calculate jacobian for estimating the variance of pij parameters
    jpij <- array( 0 , dim = c( nrow(bigNij) , ncol(bigNij) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      jpij[i,j] <-
        - bigNij[i,j] / pij[i,j]^2 -
        bigRi[i] * ( 1 - rho[j] )^2 / rowSums( rhocpij )[i]^2 -
        bigCj[j] * eta[i]^2 / colSums( nipij )[j]^2
      bigM * ( tau[j] * eta[i] )^2 / sum( taunipij )^2
    }

    # divide u_pij by the jacobian
    u_pij <- sweep( u_pij , 2:3 , jpij , "/" )

    ### muij

    # calculate linearized variable
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    ### gamma

    # calculate linearized variable
    u_gamma <- sapply( seq_len( ncol( pij ) ) , function(j) {
      rowSums( sweep( u_pij[,,j] , 2 , eta , "*" ) + sweep( u_eta , 2 , pij[,j] , "*" ) )
    } )

    ### combine results

    # return linearized variables
    reslin <-
      list(
        "psi" = u_psi ,
        "rho" = u_rho ,
        "tau" = u_tau ,
        "eta" = u_eta ,
        "pij" = u_pij ,
        "muij" = u_muij ,
        "gamma" = u_gamma )

  } )

  ### calculate variances

  # response/nonresponse probabilities
  if ( rp.variance & res$model != "MCAR" ) {

    ### psi

    # calculate variance of psi
    psi_var <- survey::svyrecvar( sweep( as.matrix( u_psi ) , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    # psi_var <- diag( psi_var )

    ### rho

    # calculate variance of rho
    rho_var <- survey::svyrecvar( sweep( as.matrix( u_rho ) , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    # rho_var <- diag( rho_var )

    ### tau

    # calculate variance of tau
    tau_var <- survey::svyrecvar( sweep( as.matrix( u_tau ) , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
    # tau_var <- diag( tau_var )

  } else {

    psi_var <- as.matrix( NA , length(psi))
    rho_var <- as.matrix( NA , length(rho))
    tau_var <- as.matrix( NA , length(tau))

  }

  ### eta

  # calculate variance of eta
  eta_var <- survey::svyrecvar( sweep( u_eta , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  # eta_var <- diag( eta_var )

  ### pij

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

  # create matrix of linearized variables
  lin_muij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_muij[,z,] ) )

  # adjust for structural zeros in transition matrix
  if ( !is.null( pij.zero ) ) {
    lin_muij[ , c( t( pij.zero ) ) ] <- 0
  }

  # calculate variance of muij
  muij_var <- survey::svyrecvar( ww * lin_muij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  muij_var <- matrix( diag( muij_var ) , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

  ### gamma

  # calculate variance of gamma
  gamma_var <- survey::svyrecvar( sweep( u_gamma , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  # gamma_var <- diag( gamma_var )

  ### net flow

  # variance
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

  if ( res$model != "MCAR" ) {

    ### calculate rao-scott correction

    # combine objects
    u_R <- yy[,,1] * ( 1 - zz[,2] )
    u_C <- yy[,,2] * ( 1 - zz[,1] )
    u_M <- ( 1 - zz[,1] ) * ( 1 - zz[,2] )
    lmat <- y1y2
    lmat <- abind::abind( list( lmat , u_R ) , along = 3 )
    lmat <- abind::abind( list( lmat , cbind( u_C , u_M , deparse.level = 0 ) ) , along = 2 )

    # variance under SRS
    smalln <- sum( ww > 0 )
    mean2 <- colSums( ww * do.call( cbind , lapply( seq_len( ncol( bigNij )+1 ) , function( z ) lmat[,z,] ) ) ) / N
    Dmat <- diag(mean2)
    iDmat <- diag(ifelse(mean2 == 0, 0, 1/mean2))
    Vsrs <- (Dmat - outer(mean2, mean2))/(smalln-1)

    # apply ratio linearization
    lmat <- apply( lmat , 2:3 , function(z) z/N - sum( ww * z ) / N^2 )

    # linearized matrix
    lmat <- do.call( cbind , lapply( seq_len( ncol( bigNij )+1 ) , function( z ) lmat[,z,] ) )

    # calculate variance
    Vdes <- survey::svyrecvar( sweep( lmat , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    # delta matrix
    Delta <- MASS::ginv( Vsrs ) %*% Vdes

    # eigenvalues of delta matrix
    smalld <- eigen( Delta , only.values = TRUE )$values[ - NCOL(Delta) ]
    smalld <- Re( smalld )

    # # first-order correction
    # meand <- mean( smalld )
    # adj.chi <- chiscore / meand

    # second-order correction
    meand <- mean( smalld )
    smalla2 <- mean( smalld^2 ) / meand^2
    pearson$statistic <- pearson$statistic/( 1 + smalla2 )
    pearson$p.value <- if ( res$model %in% c("A","B") ) pchisq( pearson$statistic , df = ( NCOL( Delta ) - 1 ) / ( 1 + smalla2 ) , lower.tail = FALSE ) else NA
    # attr( pearson$statistic, "names" ) <- "F"
    pearson$parameter <- c( df = ifelse( res$model %in% c("A","B") , ( NCOL( Delta ) - 1 ) / ( 1 + smalla2 ) , NA ) )
    pearson$method <- "Pearson's X^2: Rao & Scott adjustment"
    pearson$data.name <- "observed counts"

    # # Snedecor F test
    # meand <- mean( smalld )
    # smallc <- mean( smalld^2 ) / meand^2
    # d0 <- (ncol(Delta)/2) / smallc
    # nu <- degf( design )
    # warn <- options(warn = -1)
    # pearson$statistic <- pearson$statistic/mean(smalld)
    # pearson$p.value <- pf( pearson$statistic, d0, d0 * nu, lower.tail = FALSE )
    # attr( pearson$statistic, "names" ) <- "F"
    # pearson$parameter <- c(ndf = d0, ddf = d0 * nu)
    # pearson$method <- "Pearson's X^2: Rao & Scott adjustment"
    # pearson$data.name <- "observed counts"

  }

  # add results
  mvar[["adj.chisq"]] <- if ( res$model == "MCAR" ) NULL else pearson

  # return variance
  mvar

}
