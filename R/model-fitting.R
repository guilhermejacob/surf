# function for model fitting
ipf <- function( CountMatrix , model , tol = 1e-8 , verbose = FALSE ) {

  # Obtain sample estimates of population flows as described in table 3.1 of Rojas et al. (2014)
  bigRi <- CountMatrix[ , ncol(CountMatrix) ][ - nrow( CountMatrix ) ]
  bigCj <- CountMatrix[ nrow(CountMatrix) , ][ - ncol( CountMatrix ) ]
  bigM <- CountMatrix[ nrow( CountMatrix ) , ncol( CountMatrix ) ]
  bigNij <- as.matrix( CountMatrix[ -nrow( CountMatrix ) , -ncol( CountMatrix ) ] )
  N  <- sum( bigNij ) + sum( bigRi ) + sum( bigCj ) + bigM

  # fix zero probabilities
  if ( bigM == 0 ) {
    stop( "Models are undefined for tables without non-response in both times." )
  }

  # test if margins are not zero
  if ( any( rowSums( bigNij ) == 0 , colSums( bigNij ) == 0 ) ) stop( "a line or column total is equal to zero. consider removing the category." )

  # number of categories
  k <- ncol(bigNij)

  # number of counts
  n.counts <- (k+1)^2

  # number of model parameters
  switch( model ,
          MCAR = { n.parms <- k^2 + k } ,
          A = { n.parms <- k^2 +   k + 3 } ,
          B = { n.parms <- k^2 + 2*k + 2 } ,
          C = { n.parms <- k^2 + 3*k + 1 } ,
          D = { n.parms <- k^2 + 3*k + 1 } )

  # number of restrictions
  n.restr <- k+1

  # adjusts using onde of the models
  mfit <- switch( model , MCAR = {

    # MCAR estimation
    eta <- rowSums( bigNij ) / sum( bigNij )
    pij <- bigNij / rowSums( bigNij )

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
        "rho" = NA ,
        "tau" = NA ,
        "eta" = eta ,
        "pij" = pij ,
        "muij" = N * sweep( pij , 1 , eta , "*" ) ,
        "gamma" = colSums( sweep( pij , 1 , eta , "*" ) ) )

  } , A = {
    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # psi, rho, and tau according to Result 4.2 of Rojas (2014, p.38)
    psi <- ( sum( bigNij ) + sum( bigRi ) ) / N
    rho <- sum( bigNij ) / ( sum( bigNij ) + sum( bigRi ) )
    tau <- bigM / ( sum( bigCj ) + bigM )

    # Obtain starting values for estimating superpopulation model flow parameters
    # eta and pij according to Result 4.3 of Rojas (2014, p.45)
    eta0 <- etav <- rowSums( bigNij ) / sum( bigNij )
    pij0 <- pijv <- sweep( bigNij , 1 , rowSums( bigNij ) , "/" )
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
      if (verbose) cat( "iteration: " , stringr::str_pad( v , width = 4 , pad = " " , side = "left" ) , "\t maxdiff: " , scales::number( maxdiff , accuracy = tol ) , "\r" )

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
        "rho" = rho ,
        "tau" = tau ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) ,
        "gamma" = colSums( sweep( pijv , 1 , etav , "*" ) ) )
  } ,
  B ={
    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # rho, and tau according to Result 4.7 of Rojas (2014, p.46-47)
    rho <- sum( bigNij ) / ( sum( bigNij ) + sum( bigRi ) )
    tau <- bigM / ( sum( bigCj ) + bigM )

    # Obtain starting values for estimating superpopulation model flow parameters
    # psii, eta and pij according to Result 4.8 of Rojas (2014, p.49)
    # psii0 <- psiiv <- rep( sum( bigNij ) + sum( bigRi ) , nrow( bigNij ) ) / N
    psii0 <- psiiv <- rep( ( sum( bigNij ) + sum( bigRi ) ) / N , nrow( bigNij ) )
    eta0 <- etav <- rowSums( bigNij ) / sum( bigNij )
    pij0 <- pijv <- sweep( bigNij , 1 , rowSums( bigNij ) , "/" )
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

      # calculate changes
      these_diffs <- c( c( psiiv - psii0 ) , c( etav - eta0 ) , c( pijv - pij0 ) )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these_diffs ) )

      # process tracker
      if (verbose) cat( "iteration: " , stringr::str_pad( v , width = 4 , pad = " " , side = "left" ) , "\t maxdiff: " , scales::number( maxdiff , accuracy = tol ) , "\r" )

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
        "rho" = rho ,
        "tau" = tau ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) ,
        "gamma" = colSums( sweep( pijv , 1 , etav , "*" ) ) )

  } ,
  C = {

    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # rho, and tau according to Result 4.12 of Rojas (2014, p.54-55)
    psi <- ( sum( bigNij ) + sum( bigRi ) ) / ( sum( bigNij ) + sum( bigRi ) + sum( bigCj ) + bigM )
    rhoi <- rowSums( bigNij ) / ( rowSums( bigNij ) + bigRi )

    # Obtain starting values for estimating superpopulation model flow parameters
    # psii, eta and pij according to Result 4.13 of Rojas (2014, p.62)
    taui0 <- tauiv <- rep( bigM / ( sum( bigCj ) + bigM ) , nrow( bigNij ) )
    eta0 <- etav <- rowSums( bigNij ) / sum( bigNij )
    pij0 <- pijv <- sweep( bigNij , 1 , rowSums( bigNij ) , "/" )
    maxdiff <- Inf
    v = 0

    # Obtain maximum pseudo-likelihood estimates for superpopulation model flow parameters
    # psii, eta and pij according to Result 4.8 of Rojas (2014, p.296)
    while( maxdiff > tol ) {

      # count iteration
      v <- v+1

      # calculate auxiliary stats
      nipij <- sweep( pijv , 1 , etav , "*" )
      nirho <- etav * tauiv
      nipijrhoc <- sweep( nipij , 1 , 1 - tauiv , "*" )

      # calculate psi
      tauiv <-
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
      these_diffs <- c( c( tauiv - taui0 ) , c( etav - eta0 ) , c( pijv - pij0 ) )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these_diffs ) )

      # process tracker
      if (verbose) cat( "iteration: " , stringr::str_pad( v , width = 4 , pad = " " , side = "left" ) , "\t maxdiff: " , scales::number( maxdiff , accuracy = tol ) , "\r" )

      # test for convergence
      if ( maxdiff <= tol ) break()

      # store estimates for next iteration
      taui0 <- tauiv
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
        "rho" = rhoi ,
        "tau" = tauiv ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) ,
        "gamma" = colSums( sweep( pijv , 1 , etav , "*" ) ) )

  } ,
  D = {

    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # rho, and tau according to Result 4.16 of Rojas (2014, p.63)
    psi <- ( sum( bigNij ) + sum( bigRi ) ) / ( sum( bigNij ) + sum( bigRi ) + sum( bigCj ) +bigM )

    # Obtain starting values for estimating superpopulation model flow parameters
    # rhoj, tauj, eta and pij
    rhoj0 <- rhojv <- rep( sum( bigNij ) / ( sum( bigNij ) + sum( bigRi ) ) , ncol( bigNij ) )
    tauj0 <- taujv <- rep( bigM / ( sum( bigCj ) + bigM ) , ncol( bigNij ) )
    eta0 <- etav <- rowSums( bigNij ) / sum( bigNij )
    pij0 <- pijv <- sweep( bigNij , 1 , rowSums( bigNij ) , "/" )
    maxdiff <- Inf
    v = 0

    # Obtain maximum pseudo-likelihood estimates for superpopulation model flow parameters
    # rhoj, tauj, eta and pij according to Result 4.17 of Rojas (2014, p.64-65)
    while( maxdiff > tol ) {

      # count iteration
      v <- v+1

      # calculate auxiliary stats
      nipij <- sweep( pijv , 1 , etav , "*" )
      pijrhoc <- sweep( pijv , 2 , 1 - rhojv , "*" )
      nipijtau <- sweep( nipij , 2 , taujv , "*" )

      # calculate rhoj
      rhojv <-
        colSums( bigNij ) /
        ( colSums( bigNij ) + colSums( sweep( sweep( pijrhoc , 1 , rowSums( pijrhoc ) , "/" ) , 1 , bigRi , "*" ) ) )

      # calculate tauj
      taujv <- 1 - ( bigCj * sum( nipijtau ) ) / ( bigM * colSums( nipij ) )

      # calculate eta
      etav <-
        ( rowSums( bigNij ) + bigRi +
            rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , bigCj , "*" ) ) +
            bigM * rowSums( nipijtau ) / sum( nipijtau ) ) / N

      # calculate pij
      pijv <-
        sweep(
          bigNij +
            sweep( sweep( pijrhoc , 1 , rowSums( pijrhoc ) , "/" ) , 1 , bigRi , "*" ) +
            sweep( sweep( nipij , 2 , bigCj , "*" ) , 2 , colSums( nipij ) , "/" ) +
            bigM * ( nipijtau / sum( nipijtau ) ) , 1:2 ,
          rowSums( bigNij ) +
            rowSums( sweep( sweep( pijrhoc , 1 , rowSums( pijrhoc ) , "/" ) , 1 , bigRi , "*" ) ) +
            rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , bigCj , "*" ) ) +
            bigM * rowSums( nipijtau ) / sum( nipijtau ) , "/" )

      # calculate changes
      these_diffs <- c( c( rhojv - rhoj0 ) , c( taujv - tauj0 ) , c( etav - eta0 ) , c( pijv - pij0 ) )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these_diffs ) )

      # process tracker
      if (verbose) cat( "iteration: " , stringr::str_pad( v , width = 4 , pad = " " , side = "left" ) , "\t maxdiff: " , scales::number( maxdiff , accuracy = tol ) , "\r" )

      # test for convergence
      if ( maxdiff <= tol ) break()

      # store estimates for next iteration
      rhoj0 <- rhojv
      tauj0 <- taujv
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
        "rho" = rhojv ,
        "tau" = taujv ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) ,
        "gamma" = colSums( sweep( pijv , 1 , etav , "*" ) ) )

  } )

  # verbose treat
  if (verbose) cat( "\n" )

  ### net flows
  mfit[["delta"]] <- N * ( mfit$gamma - mfit$eta )

  # return mcar
  if ( model == "MCAR" ) return( mfit )

  ### log pseudo-likelihood
  nipij <- sweep( pijv , 1 , etav , "*" )
  psirho <- res[["psi"]] * res[["rho"]]
  psirhoc <- res[["psi"]] * ( 1 - res[["rho"]] )
  psictauc <- ( 1 - res[["psi"]] ) * ( 1 - res[["tau"]] )
  psictau <- ( 1 - res[["psi"]] ) * res[["tau"]]
  switch( model , A = {
    pll <- sum( bigNij * log( psirho * nipij ) , na.rm = TRUE ) +
      sum( bigRi * rowSums( log( psirhoc * nipij ) ) , na.rm = TRUE ) +
      sum( bigCj * colSums( log( psictauc * nipij ) ) , na.rm = TRUE ) +
      bigM * log( sum( psictau * nipij , na.rm = TRUE ) )
  } , B = {
    pll <- sum( bigNij * log( sweep( nipij , 1 , psirho , "*" ) ) , na.rm = TRUE ) +
      sum( bigRi * log( rowSums( sweep( nipij , 1 , psirhoc , "*" ) ) ) , na.rm = TRUE ) +
      sum( bigCj * log( colSums( sweep( nipij , 1 , psictauc , "*" ) ) ) , na.rm = TRUE ) +
      bigM * log( sum( sweep( nipij , 1 , psictau , "*" ) , na.rm = TRUE ) )
  } , C = {
    pll <- sum( bigNij * log( sweep( nipij , 1 , psirho , "*" ) ) , na.rm = TRUE ) +
      sum( bigRi * log( rowSums( sweep( nipij , 1 , psirhoc , "*" ) ) ) , na.rm = TRUE ) +
      sum( bigCj * log( colSums( sweep( nipij , 1 , psictauc , "*" ) ) ) , na.rm = TRUE ) +
      bigM * log( sum( sweep( nipij , 1 , psictau , "*" ) , na.rm = TRUE ) )
  } , D = {
    pll <- sum( bigNij * log( sweep( nipij , 2 , psirho , "*" ) ) , na.rm = TRUE ) +
      sum( bigRi * log( rowSums( sweep( nipij , 2 , psirhoc , "*" ) ) ) , na.rm = TRUE ) +
      sum( bigCj * log( colSums( sweep( nipij , 2 , psictauc , "*" ) ) ) , na.rm = TRUE ) +
      bigM * log( sum( sweep( nipij , 2 , psictau , "*" ) , na.rm = TRUE ) )
  } )
  mfit[["ll"]] <- pll / N

  ### unadjusted chi-distances
  observed.props  <- cbind( rbind( mfit$bigNij , mfit$bigRi ) , c( mfit$bigCj , mfit$bigM ) ) / mfit$N
  estimated.Nij   <- sweep( mfit$pij , 1 , mfit$rho * mfit$psi * mfit$eta , "*" )
  estimated.Ri    <- mfit$psi * ( 1 - mfit$rho ) * rowSums( sweep( mfit$pij , 1 , mfit$eta , "*" ) )
  estimated.Cj    <- ( 1 - mfit$psi ) * ( 1 - mfit$tau ) * colSums( sweep( mfit$pij , 1 , mfit$eta , "*" ) )
  estimated.M     <- sum( sweep( mfit$pij , 1 , mfit$tau * ( 1 - mfit$psi ) * mfit$eta , "*" ) )
  estimated.props <- cbind( rbind( estimated.Nij , estimated.Ri ) , c( estimated.Cj , estimated.M ) )
  chimat <- ( observed.props - estimated.props )^2 / estimated.props
  dimnames( chimat ) <- NULL

  # store estimated counts
  mfit[["estimated.counts"]] <- estimated.props * N

  # calculate unadjusted test score
  chiscore <- sum( N * chimat )

  # store unadjusted chi-square test score
  warn <- options(warn = -1)
  mfit[["unadj.chisq"]] <- chisq.test( matrix( 10, ncol = 3 , nrow = 3 ) , correct = FALSE )
  mfit[["unadj.chisq"]]$statistic[[1]] <- chiscore

  # recalculate p-value
  if (model %in% c( "A" , "B" ) ) {
    mfit[["unadj.chisq"]]$parameter[[1]] <- n.counts + n.restr - n.parms
    mfit[["unadj.chisq"]]$p.value <- pchisq( mfit[["unadj.chisq"]]$statistic , n.counts + n.restr - n.parms , lower.tail = FALSE )
  } else {
    mfit[["unadj.chisq"]]$parameter[[1]] <- 0
    mfit[["unadj.chisq"]]$p.value <-  NA
  }
  mfit[["unadj.chisq"]]$data.name <- "observed vs. expected counts"
  mfit[["unadj.chisq"]]$method <- "Unadjusted Pearson's Chi-squared test"

  # store countss
  mfit[["model.info"]] <- c( n.counts , n.restr , n.parms )
  mfit[["observed.counts"]] <- CountMatrix

  # return fit
  return( mfit )

}

