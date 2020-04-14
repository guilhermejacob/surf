# function for model fitting
ipf <- function( xx , ww , model , tol = 1e-8 , verbose = FALSE , starting.values = list( psi = NULL , rhoRR = NULL , rhoMM = NULL , eta = NULL , pij = NULL ) ) {

  # Obtain sample estimates of population flows as described in table 3.1 of Rojas et al. (2014)
  bigNij <- stats::xtabs( c(ww,0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE , drop.unused.levels = FALSE )
  bigRi <- bigNij[ , ncol(bigNij) ][ - nrow( bigNij ) ]
  bigCj <- bigNij[ nrow(bigNij) , ][ - ncol( bigNij ) ]
  bigM <- bigNij[ nrow( bigNij ) , ncol( bigNij ) ]
  bigNij <- as.matrix( bigNij[ -nrow( bigNij ) , -ncol( bigNij ) ] )
  N  <- sum( bigNij ) + sum( bigRi ) + sum( bigCj ) + bigM

  # ajusts using onde of the models
  mfit <- switch( model , A = {
    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # psi, rhoRR, and rhoMM according to Result 4.2 of Rojas (2014, p.38)
    psi <- ( sum( bigNij ) + sum( bigRi ) ) / N
    rhoRR <- sum( bigNij ) / ( sum( bigNij ) + sum( bigRi ) )
    rhoMM <- bigM / ( sum( bigCj ) + bigM )

    # Obtain starting values for estimating superpopulation model flow parameters
    # eta and pij according to Result 4.3 of Rojas (2014, p.45)
    eta0 <- etav <- if ( is.null( starting.values[["eta"]] ) ) rowSums( bigNij ) / sum( bigNij ) else starting.values[["eta"]]
    pij0 <- pijv <- if ( is.null( starting.values[["pij"]] ) ) sweep( bigNij , 1 , rowSums( bigNij ) , "/" ) else starting.values[["pij"]]
    maxdiff <- Inf
    v = 0

    # Obtain maximum pseudo-likelihood estimates for superpopulation model flow parameters
    # eta and pij according to Result 4.3 of Rojas (2014, p.41)
    while( maxdiff > tol ) {

      # count iteration
      v <- v+1

      # calculate auxiliary stats
      nipij <- sweep( pijv , 1 , etav , "*" )

      # calculate eta v+1
      etav <-
        ( rowSums( bigNij ) + bigRi + rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , bigCj , "*" ) ) ) /
        ( sum( bigNij ) + sum( bigRi ) + sum( bigCj ) )

      # calculate pij v+1
      pijv <-
        sweep( bigNij + sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , bigCj , "*" ) , 1 ,
               rowSums( bigNij ) + rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , bigCj , "*" ) ) , "/" )

      # calculate maximum absolute difference
      maxdiff <- max( abs( c( c( etav - eta0 ) , c( pijv - pij0 ) ) ) )

      # process tracker
      if (verbose) cat( ifelse(v==1,"\n","") , " iteration: " , stringr::str_pad( v , width = 4 , pad = " " , side = "left" ) , "\t maxdiff: " , scales::number( maxdiff , accuracy = tol ) , "\r" )

      # test for convergence
      if ( maxdiff <= tol ) break()

      # store estimates for next iteration
      eta0 <- etav
      pij0 <- pijv

    }

    # return final estimates
    res <-
      list(
        "model" = model ,
        "iter" = v ,
        "N" = N ,
        "bigNij" = bigNij ,
        "bigRi" = bigRi ,
        "bigCj" = bigCj ,
        "bigM" = bigM ,
        "psi" = psi ,
        "rhoRR" = rhoRR ,
        "rhoMM" = rhoMM ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) )

    # return list of results
    return( res )

  } ,
  B ={
    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # rhoRR, and rhoMM according to Result 4.7 of Rojas (2014, p.46-47)
    rhoRR <- sum( bigNij ) / ( sum( bigNij ) + sum( bigRi ) )
    rhoMM <- bigM / ( sum( bigCj ) + bigM )

    # Obtain starting values for estimating superpopulation model flow parameters
    # psii, eta and pij according to Result 4.8 of Rojas (2014, p.49)
    # psii0 <- psiiv <- rep( sum( bigNij ) + sum( bigRi ) , nrow( bigNij ) ) / N
    psii0 <- psiiv <- if ( is.null( starting.values[["psi"]] ) ) ( rowSums( bigNij ) + sum( bigRi ) ) / N else starting.values[["psi"]]
    eta0 <- etav <- if ( is.null( starting.values[["eta"]] ) ) rowSums( bigNij ) / sum( bigNij ) else starting.values[["eta"]]
    pij0 <- pijv <- if ( is.null( starting.values[["pij"]] ) ) sweep( bigNij , 1 , rowSums( bigNij ) , "/" ) else starting.values[["pij"]]
    maxdiff <- Inf
    v = 0

    # Obtain maximum pseudo-likelihood estimates for superpopulation model flow parameters
    # psii, eta and pij according to Result 4.8 of Rojas (2014, p.296)
    while( maxdiff > tol ) {

      # count iteration
      v <- v+1

      # calculate auxiliary stats
      nipij <- sweep( pijv , 1 , etav , "*" )
      psicni <- ( 1 - psiiv ) * etav
      psicnipij <- sweep( pijv , 1 , psicni , "*" )

      # calculate psi
      psiiv <-
        ( rowSums( bigNij ) + bigRi ) /
        ( rowSums( bigNij ) + bigRi +
            rowSums( sweep( sweep( psicnipij , 2 , colSums( psicnipij ) , "/" ) , 2 , bigCj , "*" ) ) +
            bigM * ( psicni / sum( psicni ) ) )

      # calculate eta v+1
      etav <-
        ( rowSums( bigNij ) + bigRi +
            rowSums( sweep( sweep( psicnipij , 2 , colSums( psicnipij ) , "/" ) , 2 , bigCj , "*" ) ) +
            bigM * ( psicni / sum( psicni ) ) ) / ( sum( bigNij ) + sum( bigRi ) + sum( bigCj ) + bigM )

      # calculate pij v+1
      pijv <-
        sweep( bigNij + sweep( sweep( psicnipij , 2 , colSums( psicnipij ) , "/" ) , 2 , bigCj , "*" ) , 1 ,
               rowSums( bigNij ) + rowSums( sweep( sweep( psicnipij , 2 , colSums( psicnipij ) , "/" ) , 2 , bigCj , "*" ) ) , "/" )

      # calculate changes
      these_diffs <- c( c( psiiv - psii0 ) , c( etav - eta0 ) , c( pijv - pij0 ) )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these_diffs ) )

      # process tracker
      if (verbose) cat( ifelse(v==1,"\n","") , " iteration: " , stringr::str_pad( v , width = 4 , pad = " " , side = "left" ) , "\t maxdiff: " , scales::number( maxdiff , accuracy = tol ) , "\r" )

      # test for convergence
      if ( maxdiff <= tol ) break()

      # store estimates for next iteration
      psii0 <- psiiv
      eta0 <- etav
      pij0 <- pijv

    }

    # return final estimates
    res <-
      list(
        "model" = model ,
        "iter" = v ,
        "N" = N ,
        "bigNij" = bigNij ,
        "bigRi" = bigRi ,
        "bigCj" = bigCj ,
        "bigM" = bigM ,
        "psi" = psiiv ,
        "rhoRR" = rhoRR ,
        "rhoMM" = rhoMM ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) )

    # return list of results
    return( res )
  } ,
  C = {

    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # rhoRR, and rhoMM according to Result 4.12 of Rojas (2014, p.54-55)
    psi <- ( sum( bigNij ) + sum( bigRi ) ) / ( sum( bigNij ) + sum( bigRi ) + sum( bigCj ) + bigM )
    rhoRRi <- rowSums( bigNij ) / ( rowSums( bigNij ) + bigRi )

    # Obtain starting values for estimating superpopulation model flow parameters
    # psii, eta and pij according to Result 4.13 of Rojas (2014, p.62)
    rhoMMi0 <- rhoMMiv <- if ( is.null( starting.values[["rhoMM"]] ) ) rep( bigM / ( sum( bigCj ) + bigM ) , nrow( bigNij ) ) else starting.values[["rhoMM"]]
    eta0 <- etav <- if ( is.null( starting.values[["eta"]] ) ) rowSums( bigNij ) / sum( bigNij ) else starting.values[["eta"]]
    pij0 <- pijv <- if ( is.null( starting.values[["pij"]] ) ) sweep( bigNij , 1 , rowSums( bigNij ) , "/" ) else starting.values[["pij"]]
    maxdiff <- Inf
    v = 0

    # Obtain maximum pseudo-likelihood estimates for superpopulation model flow parameters
    # psii, eta and pij according to Result 4.8 of Rojas (2014, p.296)
    while( maxdiff > tol ) {

      # count iteration
      v <- v+1

      # calculate auxiliary stats
      nipij <- sweep( pijv , 1 , etav , "*" )
      nirho <- etav * rhoMMiv
      nipijrhoc <- sweep( nipij , 1 , 1 - rhoMMiv , "*" )

      # calculate psi
      rhoMMiv <-
        ( bigM * nirho / sum( nirho ) ) /
        ( rowSums( sweep( sweep( nipijrhoc , 2 , colSums( nipijrhoc ) , "/" ) , 2 , bigCj , "*" ) ) + bigM * nirho / sum( nirho ) )

      # calculate eta v+1
      etav <-
        ( rowSums( bigNij ) + bigRi +
            rowSums( sweep( sweep( nipijrhoc , 2 , colSums( nipijrhoc ) , "/" ) , 2 , bigCj , "*" ) ) +
            bigM * ( nirho / sum( nirho ) ) ) / N

      # calculate pij v+1
      pijv <- sweep(
        bigNij + sweep( sweep( nipijrhoc , 2 , colSums( nipijrhoc ) , "/" ) , 2 , bigCj , "*" ) , 1 ,
        rowSums( bigNij ) + rowSums( sweep( sweep( nipijrhoc , 2 , colSums( nipijrhoc ) , "/" ) , 2 , bigCj , "*" ) ) , "/" )

      # calculate changes
      these_diffs <- c( c( rhoMMiv - rhoMMi0 ) , c( etav - eta0 ) , c( pijv - pij0 ) )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these_diffs ) )

      # process tracker
      if (verbose) cat( ifelse(v==1,"\n","") , " iteration: " , stringr::str_pad( v , width = 4 , pad = " " , side = "left" ) , "\t maxdiff: " , scales::number( maxdiff , accuracy = tol ) , "\r" )

      # test for convergence
      if ( maxdiff <= tol ) break()

      # store estimates for next iteration
      rhoMMi0 <- rhoMMiv
      eta0 <- etav
      pij0 <- pijv

    }

    # return final estimates
    res <-
      list(
        "model" = model ,
        "iter" = v ,
        "N" = N ,
        "bigNij" = bigNij ,
        "bigRi" = bigRi ,
        "bigCj" = bigCj ,
        "bigM" = bigM ,
        "psi" = psi ,
        "rhoRR" = rhoRRi ,
        "rhoMM" = rhoMMiv ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) )

    # return list of results
    return( res )
  } ,
  D = {

    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # rhoRR, and rhoMM according to Result 4.16 of Rojas (2014, p.63)
    psi <- ( sum( bigNij ) + sum( bigRi ) ) / ( sum( bigNij ) + sum( bigRi ) + sum( bigCj ) +bigM )

    # Obtain starting values for estimating superpopulation model flow parameters
    # rhoRRj, RhoMMj, eta and pij
    rhoRRj0 <- rhoRRjv <- if ( is.null( starting.values[["rhoRR"]] ) ) colSums( bigNij ) / ( colSums( bigNij ) + sum( bigRi ) ) else starting.values[["rhoRR"]]
    rhoMMj0 <- rhoMMjv <- if ( is.null( starting.values[["rhoMM"]] ) ) 1 - bigCj / bigM else starting.values[["rhoMM"]]
    eta0 <- etav <- if ( is.null( starting.values[["eta"]] ) ) rowSums( bigNij ) / sum( bigNij ) else starting.values[["eta"]]
    pij0 <- pijv <- if ( is.null( starting.values[["pij"]] ) ) sweep( bigNij , 1 , rowSums( bigNij ) , "/" ) else starting.values[["pij"]]
    maxdiff <- Inf
    v = 0

    # Obtain maximum pseudo-likelihood estimates for superpopulation model flow parameters
    # rhoRRj, RhoMMj, eta and pij according to Result 4.17 of Rojas (2014, p.64-65)
    while( maxdiff > tol ) {

      # count iteration
      v <- v+1

      # calculate auxiliary stats
      nipij <- sweep( pijv , 1 , etav , "*" )
      pijrhoRRc <- sweep( pijv , 2 , 1 - rhoRRjv , "*" )
      nipijrhoMM <- sweep( nipij , 2 , rhoMMjv , "*" )

      # calculate rhoRRj
      rhoRRjv <-
        colSums( bigNij ) /
        ( colSums( bigNij ) + colSums( sweep( sweep( pijrhoRRc , 1 , rowSums( pijrhoRRc ) , "/" ) , 1 , bigRi , "*" ) ) )

      # calculate rhoMMj
      rhoMMjv <- 1 - ( bigCj * sum( nipijrhoMM ) ) / ( bigM * colSums( nipij ) )

      # calculate eta
      etav <-
        ( rowSums( bigNij ) + bigRi +
            rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , bigCj , "*" ) ) +
            bigM * rowSums( nipijrhoMM ) / sum( nipijrhoMM ) ) / N

      # calculate pij
      pijv <-
        sweep(
          bigNij +
            sweep( sweep( pijrhoRRc , 1 , rowSums( pijrhoRRc ) , "/" ) , 1 , bigRi , "*" ) +
            sweep( sweep( nipij , 2 , bigCj , "*" ) , 2 , colSums( nipij ) , "/" ) +
            bigM * ( nipijrhoMM / sum( nipijrhoMM ) ) , 1:2 ,
          rowSums( bigNij ) +
            rowSums( sweep( sweep( pijrhoRRc , 1 , rowSums( pijrhoRRc ) , "/" ) , 1 , bigRi , "*" ) ) +
            rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , bigCj , "*" ) ) +
            bigM * rowSums( nipijrhoMM ) / sum( nipijrhoMM ) , "/" )

      # calculate changes
      these_diffs <- c( c( rhoRRjv - rhoRRj0 ) , c( rhoMMjv - rhoMMj0 ) , c( etav - eta0 ) , c( pijv - pij0 ) )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these_diffs ) )

      # process tracker
      if (verbose) cat( ifelse(v==1,"\n","") , " iteration: " , stringr::str_pad( v , width = 4 , pad = " " , side = "left" ) , "\t maxdiff: " , scales::number( maxdiff , accuracy = tol ) , "\r" )

      # test for convergence
      if ( maxdiff <= tol ) break()

      # store estimates for next iteration
      rhoRRj0 <- rhoRRjv
      rhoMMj0 <- rhoMMjv
      eta0 <- etav
      pij0 <- pijv

    }

    # return final estimates
    res <-
      list(
        "model" = model ,
        "iter" = v ,
        "N" = N ,
        "bigNij" = bigNij ,
        "bigRi" = bigRi ,
        "bigCj" = bigCj ,
        "bigM" = bigM ,
        "psi" = psi ,
        "rhoRR" = rhoRRjv ,
        "rhoMM" = rhoMMjv ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) )

    # return list of results
    return( res )
  } )

  # verbose treat
  if (verbose) cat( "\n\n" )

  # return fit
  mfit

}

# function for model variance calculation
ipf_variance <- function( xx , ww , res , model , design ) {

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
  mvar <- switch( res$model , A = {

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
    a4 <- sum( bigNij ) / ( sum( bigNij ) + sum( bigRi ) )^2

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

    # calculate variance of u_pij sum
    lin_pij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_pij[,z,] ) )
    pij_var <- diag( survey::svyrecvar( ww * lin_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    pij_var <- matrix( pij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    # calculate variance of muij
    lin_muij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_muij[,z,] ) )
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
    a4 <- sum( bigNij ) / ( sum( bigNij ) + sum( bigRi ) )^2

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
        ( 1- zz[,1] ) * ( 1 - zz[,2] ) * rowSums( nipij )[i] / sum( psicnipij ) # wrong index in the text
    }

    # Calculate jacobian for estimating the variance of psi parameters
    # not sure about 5.77
    jpsi <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jpsi[i] <-
        - sum( ww * rowSums( y1y2[,i,] ) / psi[i]^2 ) -
        # sum( ww  * ( 1 - zz[,1] ) * rowSums( sweep( yy[,,2] , 2 , ( nipij[i,] / colSums( psicnipij )^2 ) , "*" ) ) ) - # 1st attempt
        sum( ww  * ( 1 - zz[,1] ) * rowSums( sweep( yy[,,2] , 2 , ( nipij[i,] / colSums( sweep( nipij , 1 , (1 - psi)^2 , "*" ) ) ) , "*" ) ) ) - # 2nd attempt
        sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * rowSums( nipij / sum( psicnipij )^2 )[i] )
      # sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * rowSums( nipij )[i] / sum( sweep( nipij , 1 , (1 - psi)^2 , "*" ) ) )
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
        ( rowSums( y1y2[,i,] ) + yy[,i,1] * ( 1 - zz[,2] ) ) / eta[i] +
        rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , psicpij[i,] / colSums( psicnipij ) , "*" ) ) +
        ( 1- zz[,1] ) * ( 1 - zz[,2] ) * rowSums( psicpij )[i] / sum( psicnipij ) - 1
    }

    # Calculate jacobian for estimating the variance of eta parameters
    jeta <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jeta[i] <-
        - sum( ww * y1y2[,i,] ) / eta[i]^2 -
        sum( ww * yy[,i,1] * ( 1 - zz[,2] ) ) / eta[i]^2 -
        sum( ww * ( 1 - zz[,1] ) * rowSums( sweep( yy[,,2] , 2 , psicpij[i,]^2 / colSums( psicnipij )^2 , "*" ) ) ) -
        sum( ww * ( 1- zz[,1] ) * ( 1 - zz[,2] ) * rowSums( psicpij )[i]^2 / sum( psicnipij )^2 )
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

    # calculate variance of u_pij sum
    lin_pij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_pij[,z,] ) )
    pij_var <- diag( survey::svyrecvar( ww * lin_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    pij_var <- matrix( pij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    # calculate variance of muij
    lin_muij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_muij[,z,] ) )
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

    # special variables
    a3 <- sum( bigRi ) / ( sum( bigNij ) + sum( bigRi ) )^2
    a4 <- sum( bigNij ) / ( sum( bigNij ) + sum( bigRi ) )^2

    # linearized variable e_psi
    e_rhoRR <- a3 + a4 * ( 1 - zz[,2] )

    # calculate variance
    rhoRR_var <- survey::svyrecvar( ww * e_rhoRR , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

    ### rhoMM

    # Calculate scores for estimating the variance of rhoMM parameters
    u_rhoMM <- array( 0 , dim = c( nrow(xx) , nrow( bigNij ) ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      u_rhoMM[,i] <-
        - rowSums( sweep( yy[,,2] * ( 1 - zz[,1] ) , 2 , nipij[i,] / colSums( rhocnipij ) ) ) +
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * sum( nipij[i,] ) / sum( rhonipij )
    }

    # Calculate jacobian for estimating the variance of rhoMM parameters
    jrhoMM <- vector( "numeric" , length = nrow( bigNij ) )
    for ( i in seq_len( nrow( bigNij ) ) ) {
      jrhoMM[i] <-
        - sum( ww * ( 1 - zz[,1] ) * sweep( yy[,,2] , 2 , nipij[i,]^2 / colSums( rhocnipij )^2 , "*" ) ) - sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * sum( nipij[i,] )^2 ) / sum( rhonipij )^2
      # + sum( ww * ( 1 - zz[,1] ) * sweep( yy[,,2] , 2 , nipij[i,] / colSums( sweep( nipij , 1 , ( 1 - rhoMM )^2 , "*" ) ) , "*" ) ) - sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * sum( nipij[i,] ) ) / sum( sweep( nipij , 1 , rhoMM^2 , "*" ) ) #
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

    # calculate variance of u_pij sum
    lin_pij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_pij[,z,] ) )
    pij_var <- diag( survey::svyrecvar( ww * lin_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    pij_var <- matrix( pij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    # calculate variance of muij
    lin_muij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_muij[,z,] ) )
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
        rowSums( y1y2[,,j] / rhoRR[j] ) -
        rowSums( sweep( yy[,,1] * ( 1 - zz[,2] ) , 2 , pij[,j] / rowSums( rhoRRcpij ) , "*" ) )
    }

    # Calculate jacobian for estimating the variance of rhoRR parameters
    jrhoRR <- vector( "numeric" , length = nrow( bigNij ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      jrhoRR[j] <-
        - sum( ww * rowSums( y1y2[,,j] ) ) / rhoRR[j]^2 -
        sum( ww * ( 1 - zz[,2] ) * rowSums( sweep( yy[,,1] , 2 , pij[,j]^2 / rowSums( rhoRRcpij )^2 , "*" ) ) )
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
        ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * colSums( nipij )[j] / sum( rhoMMnipij )
    }

    # Calculate jacobian for estimating the variance of rhoMM parameters
    jrhoMM <- vector( "numeric" , length = nrow( bigNij ) )
    for ( j in seq_len( nrow( bigNij ) ) ) {
      jrhoMM[j] <-
        sum( ww * yy[,j,2] * ( 1 - zz[,1] ) ) / ( 1 - rhoMM[j] )^2 -
        sum( ww * ( 1 - zz[,1] ) * ( 1 - zz[,2] ) ) * colSums( nipij )[j]^2 / sum( rhoMMnipij )^2
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

    # calculate variance of u_pij sum
    lin_pij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_pij[,z,] ) )
    pij_var <- diag( survey::svyrecvar( ww * lin_pij , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )
    pij_var <- matrix( pij_var , nrow = nrow( bigNij ) , ncol = ncol( bigNij ) , byrow = TRUE )

    ### muij

    # calculate variance
    u_muij <- array( 0 , dim = c( nrow( xx ) , nrow( bigNij ) , ncol( bigNij ) ) )
    for ( i in seq_len( nrow(bigNij) ) ) for ( j in seq_len( ncol( bigNij ) ) ) {
      u_muij[,i,j] <- nipij[i,j] + N * ( pij[i,j] * u_eta[,i] + eta[i] * u_pij[,i,j] )
    }

    # calculate variance of muij
    lin_muij <- do.call( cbind , lapply( seq_len(ncol( bigNij)) , function( z ) u_muij[,z,] ) )
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

