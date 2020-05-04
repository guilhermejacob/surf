# function for model fitting
ipf <- function( xx , ww , model , tol = 1e-8 , verbose = FALSE , starting.values = list( "psi" = NULL , "rhoRR" = NULL , "rhoMM" = NULL , "eta" = NULL , "pij" = NULL ) ,  pij.zero = NULL ) {

  # Obtain sample estimates of population flows as described in table 3.1 of Rojas et al. (2014)
  bigNij <- stats::xtabs( c(ww,0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE , drop.unused.levels = FALSE )
  bigRi <- bigNij[ , ncol(bigNij) ][ - nrow( bigNij ) ]
  bigCj <- bigNij[ nrow(bigNij) , ][ - ncol( bigNij ) ]
  bigM <- bigNij[ nrow( bigNij ) , ncol( bigNij ) ]
  bigNij <- as.matrix( bigNij[ -nrow( bigNij ) , -ncol( bigNij ) ] )
  N  <- sum( bigNij ) + sum( bigRi ) + sum( bigCj ) + bigM

  # test if margins are not zero
  if ( any( rowSums( bigNij ) == 0 , colSums( bigNij ) == 0 ) ) stop( "a line or column total is equal to zero. consider removing the category." )

  # ajusts using onde of the models
  mfit <- switch( model , MCAR = {

    # MCAR estimation
    eta <- rowSums( bigNij ) / sum( bigNij )
    pij <- bigNij / rowSums( bigNij )

    # adjust for structural zeros in transition matrix
    if ( !is.null( pij.zero ) ) pij[ pij.zero ] <- 0

    # return final estimates
    res <-
      list(
        "model" = model ,
        "iter" = NA ,
        "N" = N ,
        "bigNij" = bigNij ,
        "bigRi" = bigRi ,
        "bigCj" = bigCj ,
        "bigM" = bigM ,
        "psi" = NA ,
        "rhoRR" = NA ,
        "rhoMM" = NA ,
        "eta" = eta ,
        "pij" = pij ,
        "muij" = N * sweep( pij , 1 , eta , "*" ) ,
        "pij.zero" = pij.zero )

    # return list of results
    return( res )

  } , A = {
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

      # adjust for structural zeros in transition matrix
      if ( !is.null( pij.zero ) ) pij0[ pij.zero ] <- pijv[ pij.zero ] <- 0

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
        "muij" = N * sweep( pijv , 1 , etav , "*" ) ,
        "pij.zero" = pij.zero )

  } ,
  B ={
    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # rhoRR, and rhoMM according to Result 4.7 of Rojas (2014, p.46-47)
    rhoRR <- sum( bigNij ) / ( sum( bigNij ) + sum( bigRi ) )
    rhoMM <- bigM / ( sum( bigCj ) + bigM )

    # Obtain starting values for estimating superpopulation model flow parameters
    # psii, eta and pij according to Result 4.8 of Rojas (2014, p.49)
    # psii0 <- psiiv <- rep( sum( bigNij ) + sum( bigRi ) , nrow( bigNij ) ) / N
    psii0 <- psiiv <- if ( is.null( starting.values[["psi"]] ) ) rep( ( sum( bigNij ) + sum( bigRi ) ) / N , nrow( bigNij ) ) else starting.values[["psi"]]
    eta0 <- etav <- if ( is.null( starting.values[["eta"]] ) ) rowSums( bigNij ) / sum( bigNij ) else starting.values[["eta"]]
    pij0 <- pijv <- if ( is.null( starting.values[["pij"]] ) ) sweep( bigNij , 2 , rowSums( bigNij ) , "/" ) else starting.values[["pij"]]
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
      psicpij <- sweep( pijv , 1 , 1 - psiiv , "*" )

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
            bigM * ( psicni / sum( psicni ) ) ) / N

      # calculate pij v+1
      pijv <-
        sweep( bigNij + sweep( sweep( psicnipij , 2 , colSums( psicnipij ) , "/" ) , 2 , bigCj , "*" ) , 1 ,
               rowSums( bigNij ) + rowSums( sweep( sweep( psicnipij , 2 , colSums( psicnipij ) , "/" ) , 2 , bigCj , "*" ) ) , "/" )

      # adjust for structural zeros in transition matrix
      if ( !is.null( pij.zero ) ) pij0[ pij.zero ] <- pijv[ pij.zero ] <- 0

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
        "muij" = N * sweep( pijv , 1 , etav , "*" ) ,
        "pij.zero" = pij.zero )

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

      # adjust for structural zeros in transition matrix
      if ( !is.null( pij.zero ) ) pij0[ pij.zero ] <- pijv[ pij.zero ] <- 0

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
        "muij" = N * sweep( pijv , 1 , etav , "*" ) ,
        "pij.zero" = pij.zero )

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

      # adjust for structural zeros in transition matrix
      if ( !is.null( pij.zero ) ) pij0[ pij.zero ] <- pijv[ pij.zero ] <- 0

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
        "muij" = N * sweep( pijv , 1 , etav , "*" ) ,
        "pij.zero" = pij.zero )

  } )

  ### unadjusted chi-distances
  observed.props  <- cbind( rbind( mfit$bigNij , mfit$bigRi ) , c( mfit$bigCj , mfit$bigM ) ) / mfit$N
  estimated.Nij   <- sweep( mfit$pij , 1 , mfit$rhoRR * mfit$psi * mfit$eta , "*" )
  estimated.Ri    <- mfit$psi * ( 1 - mfit$rhoRR ) * rowSums( sweep( mfit$pij , 1 , mfit$eta , "*" ) )
  estimated.Cj    <- ( 1 - mfit$psi ) * ( 1 - mfit$rhoMM ) * colSums( sweep( mfit$pij , 1 , mfit$eta , "*" ) )
  estimated.M     <- sum( sweep( mfit$pij , 1 , mfit$rhoMM * ( 1 - mfit$psi ) * mfit$eta , "*" ) )
  estimated.props <- cbind( rbind( estimated.Nij , estimated.Ri ) , c( estimated.Cj , estimated.M ) )
  chimat <- ( observed.props - estimated.props )^2 / estimated.props
  dimnames( chimat ) <- NULL

  # store estimated counts
  mfit[["estimated.counts"]] <- estimated.props * N

  # adjust for structural zeros in transition matrix
  if ( !is.null( pij.zero ) ) chimat[ pij.zero ] <- 0

  # calculate number of observations
  smalln <- sum( ww > 0 )

  # calculate unadjusted test score
  chimat <- smalln * chimat
  chiscore <- sum( chimat )

  # store unadjusted chi-square test score
  mfit[["unadj.chisq"]] <- chiscore

  # verbose treat
  if (verbose) cat( "\n\n" )

  # return fit
  mfit

}

