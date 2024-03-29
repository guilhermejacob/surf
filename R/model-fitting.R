# function for model fitting
ipf <- function( CountMatrix , model , tol = NULL , maxit = 500 , verbose = FALSE , keep.info = FALSE ) {

  # limit fitting algorithm information matrix
  if( is.infinite( maxit ) ) keep.info <- FALSE

  # Obtain sample estimates of population flows as described in table 3.1 of Rojas et al. (2014)
  Ri <- CountMatrix[ , ncol(CountMatrix) ][ - nrow( CountMatrix ) ]
  Cj <- CountMatrix[ nrow(CountMatrix) , ][ - ncol( CountMatrix ) ]
  M <- CountMatrix[ nrow( CountMatrix ) , ncol( CountMatrix ) ]
  Nij <- as.matrix( CountMatrix[ -nrow( CountMatrix ) , -ncol( CountMatrix ) ] )
  N  <- sum( Nij ) + sum( Ri ) + sum( Cj ) + M

  # set tolerance
  if ( is.null( tol ) ) tol <- 1/N

  # fix zero probabilities
  if ( M == 0 ) {
    stop( "Models are undefined for tables without non-response in both times." )
  }

  # test if margins are not zero
  if ( any( rowSums( Nij ) == 0 , colSums( Nij ) == 0 ) ) stop( "a line or column total is equal to zero. consider removing the category." )

  # number of categories
  k <- ncol(Nij)

  # number of counts
  n.counts <- (k+1)^2

  # number of model parameters
  switch( model ,
          A = { n.parms <- k^2 +   k + 3 } ,
          B = { n.parms <- k^2 + 2*k + 2 } ,
          C = { n.parms <- k^2 + 3*k + 1 } ,
          D = { n.parms <- k^2 + 3*k + 1 } )

  # number of restrictions
  n.restr <- k+1

  # define (pseudo) log-likelihood
  this.loglik <- switch( model , A = modelA.loglik , B = modelB.loglik , C = modelC.loglik , D = modelD.loglik )

  # define expected proportions function
  this.expfun <- switch( model , A = modelA.expected , B = modelB.expected , C = modelC.expected , D = modelD.expected )

  # create algorithm info
  if ( keep.info ) {
    info.mat <- matrix( as.numeric( NA ) , ncol = n.parms + 3 , nrow = maxit )
  }

  # initial settings for iterations
  maxdiff <- 1
  v <- 0
  last.loglik <- 0
  ndig <- ceiling( -log(tol,10) ) + 2

  # adjusts using one of the models
  mfit <- switch( model , A = {
    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # psi, rho, and tau according to Result 4.2 of Rojas (2014, p.38)
    psi <- ( sum( Nij ) + sum( Ri ) ) / N
    rho <- sum( Nij ) / ( sum( Nij ) + sum( Ri ) )
    tau <- M / ( sum( Cj ) + M )

    # Obtain starting values for estimating superpopulation model flow parameters
    # eta and pij according to Result 4.3 of Rojas (2014, p.45)
    eta0 <- etav <- rowSums( Nij ) / sum( Nij )
    pij0 <- pijv <- sweep( Nij , 1 , rowSums( Nij ) , "/" )

    # Obtain maximum pseudo-likelihood estimates for superpopulation model flow parameters
    # eta and pij according to Result 4.3 of Rojas (2014, p.41)
    while( maxdiff >= tol ) {

      # count iteration
      v <- v+1

      # calculate auxiliary stats
      nipij <- sweep( pij0 , 1 , eta0 , "*" )

      # calculate eta v+1
      etav <-
        ( rowSums( Nij ) + Ri + rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , Cj , "*" ) ) + M * eta0 ) / N
      # ( rowSums( Nij ) + Ri + rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , Cj , "*" ) ) ) / (N - M)

      # calculate pij v+1
      pijv <-
        sweep( Nij + sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , Cj , "*" ) , 1 ,
               rowSums( Nij ) + rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , Cj , "*" ) ) , "/" )

      # calculate differences
      these.diffs <- c( etav - eta0 , pijv - pij0 )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these.diffs ) )

      # calculate pseudo-likelihood
      iter.loglik <- this.loglik( c( psi , rho , tau , etav , t( pijv ) ) , CountMatrix )
      delta.loglik <- iter.loglik - last.loglik
      delta.loglik <- delta.loglik / N

      # process tracker
      if (verbose) cat( sprintf( paste0( "iteration: %5s\t maxdiff: %1." , ndig , "f \t loglik: %-10." , ndig , "f \n" ) , v , maxdiff , iter.loglik ) )

      # store fitting information
      if ( keep.info ) {
        info.mat[ v , ] <- c( "iter" = v , psi , rho , tau , etav , pijv , "maxdiff" = maxdiff , "pll" = iter.loglik )
      }

      # test for non-convergence
      if ( v >= maxit ) {
        stop( "Algorithm has not converged after " , v , " iterations." )
        break()
      }

      # store estimates for next iteration
      eta0 <- etav
      pij0 <- pijv
      last.loglik <- iter.loglik


    }

    # return final estimates
    res <-
      list(
        "model" = model ,
        "iter" = v ,
        "N" = N ,
        "Nij" = Nij ,
        "Ri" = Ri ,
        "Cj" = Cj ,
        "M" = M ,
        "psi" = psi ,
        "rho" = rho ,
        "tau" = tau ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) )
  } , B = {
    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # rho, and tau according to Result 4.7 of Rojas (2014, p.46-47)
    rho <- sum( Nij ) / ( sum( Nij ) + sum( Ri ) )
    tau <- M / ( sum( Cj ) + M )

    # Obtain starting values for estimating superpopulation model flow parameters
    # psi, eta and pij according to Result 4.8 of Rojas (2014, p.49)
    psi0 <- psiv <- rep( ( sum( Nij ) + sum( Ri ) ) / N , nrow( Nij ) ) # original
    # psi0 <- psiv <- ( rowSums( Nij ) + Ri ) / N # modified #1
    # psi0 <- psiv <- 1/( 1 + ( Cj + M )/( rowSums( Nij ) + Ri ) ) # modified #2
    eta0 <- etav <- rowSums( Nij ) / sum( Nij )
    pij0 <- pijv <- sweep( Nij , 1 , rowSums( Nij ) , "/" )

    # Obtain maximum pseudo-likelihood estimates for superpopulation model flow parameters
    # psi, eta and pij according to Result 4.8 of Rojas (2014, p.296)
    while( maxdiff >= tol ) {

      # count iteration
      v <- v+1

      # calculate auxiliary stats
      nipij <- sweep( pij0 , 1 , eta0 , "*" )
      psicni <- ( 1 - psi0 ) * eta0
      psicnipij <- sweep( nipij , 1 , 1 - psi0 , "*" )

      # calculate psi
      psiv <-
        ( rowSums( Nij ) + Ri ) /
        ( (rowSums( Nij ) + Ri) +
            rowSums( sweep( sweep( psicnipij , 2 , colSums( psicnipij ) , "/" ) , 2 , Cj , "*" ) ) +
            ( M * ( psicni / sum( psicni ) ) ) )

      # calculate eta v+1
      etav <-
        ( rowSums( Nij ) + Ri +
            rowSums( sweep( sweep( psicnipij , 2 , colSums( psicnipij ) , "/" ) , 2 , Cj , "*" ) ) +
            ( M * ( psicni / sum( psicni ) ) ) ) / N

      # calculate pij v+1
      pijv <- Nij + sweep( psicnipij , 2 , Cj/colSums( psicnipij ) , "*" )
      pijv <- sweep( pijv , 1 , rowSums( pijv ) , "/" )

      # calculate differences
      these.diffs <- c( c( psiv - psi0 ) , c( etav - eta0 ) , c( pijv - pij0 ) )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these.diffs ) )

      # calculate pseudo-likelihood
      iter.loglik <- this.loglik( c( psiv , rho , tau , etav , t( pijv ) ) , CountMatrix )

      # process tracker
      if (verbose) cat( sprintf( paste0( "iteration: %5s\t maxdiff: %1." , ndig , "f \t loglik: %-10." , ndig , "f \n" ) , v , maxdiff , iter.loglik ) )

      # store fitting information
      if ( keep.info ) {
        info.mat[ v , ] <- c( v , psiv , rho , tau , etav , pijv , maxdiff , iter.loglik )
      }

      # test for non-convergence
      if ( v >= maxit ) {
        stop( "Algorithm has not converged after " , v , " iterations." )
        break()
      }

      # store estimates for next iteration
      psi0 <- psiv
      eta0 <- etav
      pij0 <- pijv
      last.loglik <- iter.loglik


    }

    # return final estimates
    res <-
      list(
        "model" = model ,
        "iter" = v ,
        "N" = N ,
        "Nij" = Nij ,
        "Ri" = Ri ,
        "Cj" = Cj ,
        "M" = M ,
        "psi" = psiv ,
        "rho" = rho ,
        "tau" = tau ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) )

  } ,
  C = {

    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # rho, and tau according to Result 4.12 of Rojas (2014, p.54-55)
    psi <- ( sum( Nij ) + sum( Ri ) ) / ( sum( Nij ) + sum( Ri ) + sum( Cj ) + M )
    rhoi <- rowSums( Nij ) / ( rowSums( Nij ) + Ri )

    # Obtain starting values for estimating superpopulation model flow parameters
    # psi, eta and pij according to Result 4.13 of Rojas (2014, p.62)
    taui0 <- tauiv <- rep( M / ( sum( Cj ) + M ) , nrow( Nij ) )
    eta0 <- etav <- rowSums( Nij ) / sum( Nij )
    pij0 <- pijv <- sweep( Nij , 1 , rowSums( Nij ) , "/" )

    # Obtain maximum pseudo-likelihood estimates for superpopulation model flow parameters
    # psi, eta and pij according to Result 4.8 of Rojas (2014, p.296)
    while( maxdiff >= tol ) {

      # count iteration
      v <- v+1

      # calculate auxiliary stats
      nipij <- sweep( pijv , 1 , etav , "*" )
      nitau <- etav * tauiv
      nipijrhoc <- sweep( nipij , 1 , 1 - tauiv , "*" )

      # calculate psi
      tauiv <-
        ( M * nitau / sum( nitau ) ) /
        ( rowSums( sweep( sweep( nipijrhoc , 2 , colSums( nipijrhoc ) , "/" ) , 2 , Cj , "*" ) ) + M * nitau / sum( nitau ) )

      # calculate eta v+1
      etav <-
        ( rowSums( Nij ) + Ri +
            rowSums( sweep( sweep( nipijrhoc , 2 , colSums( nipijrhoc ) , "/" ) , 2 , Cj , "*" ) ) +
            M * ( nitau / sum( nitau ) ) ) / N

      # calculate pij v+1
      pijv <- sweep(
        Nij + sweep( sweep( nipijrhoc , 2 , colSums( nipijrhoc ) , "/" ) , 2 , Cj , "*" ) , 1 ,
        rowSums( Nij ) + rowSums( sweep( sweep( nipijrhoc , 2 , colSums( nipijrhoc ) , "/" ) , 2 , Cj , "*" ) ) , "/" )

      # check value consistency
      while ( any( tauiv < 0 | tauiv > 1 ) ) tauiv <- ( tauiv + taui0 ) / 2

      # calculate differences
      these.diffs <- c( c( tauiv - taui0 ) , c( etav - eta0 ) , c( pijv - pij0 ) )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these.diffs ) )

      # calculate pseudo-likelihood
      iter.loglik <- this.loglik( c( psi , rhoi , tauiv , etav , t( pijv ) ) , CountMatrix )

      # process tracker
      if (verbose) cat( sprintf( paste0( "iteration: %5s\t maxdiff: %1." , ndig , "f \t loglik: %-10." , ndig , "f \n" ) , v , maxdiff , iter.loglik ) )

      # store fitting information
      if ( keep.info ) {
        info.mat[ v , ] <- c( v , psi , rhoi , tauiv , etav , pijv , maxdiff , iter.loglik )
      }

      # test for non-convergence
      if ( v >= maxit ) {
        stop( "Algorithm has not converged after " , v , " iterations." )
        break()
      }

      # store estimates for next iteration
      taui0 <- tauiv
      eta0 <- etav
      pij0 <- pijv
      last.loglik <- iter.loglik


    }

    # return final estimates
    res <-
      list(
        "model" = model ,
        "iter" = v ,
        "N" = N ,
        "Nij" = Nij ,
        "Ri" = Ri ,
        "Cj" = Cj ,
        "M" = M ,
        "psi" = psi ,
        "rho" = rhoi ,
        "tau" = tauiv ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) )

  } ,
  D = {

    # issue about results
    if( any( M <= Cj ) ) {
      warning( "Cj > M: tau estimates for model D might be wrong." , call. = FALSE , immediate. = TRUE )
    }

    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # rho, and tau according to Result 4.16 of Rojas (2014, p.63)
    psi <- ( sum( Nij ) + sum( Ri ) ) / ( sum( Nij ) + sum( Ri ) + sum( Cj ) +M )

    # Obtain starting values for estimating superpopulation model flow parameters
    # rhoj, tauj, eta and pij
    rhoj0 <- rhojv <- rep( sum( Nij ) / ( sum( Nij ) + sum( Ri ) ) , ncol( Nij ) )
    # tauj0 <- taujv <- ifelse( M > Cj , 1 - Cj / M , M / ( M + Cj ) )
    tauj0 <- taujv <- M / ( M + Cj )
    # tauj0 <- taujv <- rep( M / ( sum( Cj ) + M ) , ncol( Nij ) )
    # tauj0 <- taujv <- if ( all( Cj < M ) ) ( M - Cj ) / M  else ( M - ( sum( Cj ) - Cj ) ) / ( M + Cj )
    # tauj0 <- taujv <- ifelse( M > Cj , 1 - Cj / M , M / ( M + Cj ) )
    # tauj0 <- taujv <- runif( ncol( Nij ) , .05 , .95 )
    eta0 <- etav <- rowSums( Nij ) / sum( Nij )
    pij0 <- pijv <- sweep( Nij , 1 , rowSums( Nij ) , "/" )

    # Obtain maximum pseudo-likelihood estimates for superpopulation model flow parameters
    # rhoj, tauj, eta and pij according to Result 4.17 of Rojas (2014, p.64-65)
    while( maxdiff >= tol ) {

      # count iteration
      v <- v+1

      # calculate auxiliary stats
      nipij <- sweep( pijv , 1 , etav , "*" )
      pijrhoc <- sweep( pijv , 2 , 1 - rhojv , "*" )
      nipijrhoc <- sweep( nipij , 2 , 1 - rhojv , "*" )
      nipijtau <- sweep( nipij , 2 , taujv , "*" )
      Tmat <- outer( etav , taujv )

      # calculate rhoj
      rhojv <-
        colSums( Nij ) /
        ( colSums( Nij ) + colSums( sweep( sweep( nipijrhoc , 1 , Ri , "*" ) , 1 , rowSums( nipijrhoc ) , "/" ) ) )

      # calculate tauj #1
      # taujv <- 1 - ( Cj * sum( nipijtau ) ) / ( M * colSums( nipij ) ) # Stasny
      # if ( any( taujv > 1 | taujv < 0 ) ) {
      #   Avec <- colSums( nipijtau ) / sum( nipijtau )
      #   taujv <-  (M * Avec) / ( Cj + M * Avec )
      # }
      Avec <- colSums( nipijtau ) / sum( nipijtau )
      taujv <-  (M * Avec) / ( Cj + M * Avec )

      # calculate eta
      etav <-
        ( rowSums( Nij ) + Ri +
            rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , Cj , "*" ) ) +
            M * rowSums( nipijtau ) / sum( nipijtau ) ) / N

      # calculate pij
      pijv <-
        sweep(
          Nij +
            sweep( sweep( pijrhoc , 1 , rowSums( pijrhoc ) , "/" ) , 1 , Ri , "*" ) +
            sweep( sweep( nipij , 2 , Cj , "*" ) , 2 , colSums( nipij ) , "/" ) +
            M * ( nipijtau / sum( nipijtau ) ) , 1:2 ,
          rowSums( Nij ) +
            rowSums( sweep( sweep( pijrhoc , 1 , rowSums( pijrhoc ) , "/" ) , 1 , Ri , "*" ) ) +
            rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , Cj , "*" ) ) +
            M * rowSums( nipijtau ) / sum( nipijtau ) , "/" )

      # # check value consistency #2
      # while ( any( taujv < 0 | taujv > 1 ) ) taujv <- ( taujv + tauj0 ) / 2
      # while ( any( rhojv < 0 | rhojv > 1 ) ) rhojv <- ( rhojv + rhoj0 ) / 2

      # # damping
      # taujv <- ( taujv + tauj0 ) / 2
      # rhojv <- ( rhojv + rhoj0 ) / 2

      # calculate differences
      these.diffs <- c( c( rhojv - rhoj0 ) , c( taujv - tauj0 ) , c( etav - eta0 ) , c( pijv - pij0 ) )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these.diffs ) )

      # calculate pseudo-likelihood
      iter.loglik <- this.loglik( c( psi , rhojv , taujv , etav , t( pijv ) ) , CountMatrix )

      # process tracker
      if (verbose) cat( sprintf( paste0( "iteration: %5s\t maxdiff: %1." , ndig , "f \t loglik: %-10." , ndig , "f \n" ) , v , maxdiff , iter.loglik ) )

      # store fitting information
      if ( keep.info ) {
        info.mat[ v , ] <- c( v , psi , rhojv , taujv , etav , pijv , maxdiff , iter.loglik )
      }

      # test for non-convergence
      if ( v >= maxit ) {
        stop( "Algorithm has not converged after " , v , " iterations." )
        break()
      }

      # store estimates for next iteration
      rhoj0 <- rhojv
      tauj0 <- taujv
      eta0 <- etav
      pij0 <- pijv
      last.loglik <- iter.loglik

    }

    # return final estimates
    res <-
      list(
        "model" = model ,
        "iter" = v ,
        "N" = N ,
        "Nij" = Nij ,
        "Ri" = Ri ,
        "Cj" = Cj ,
        "M" = M ,
        "psi" = psi ,
        "rho" = rhojv ,
        "tau" = taujv ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) )

  } )

  # verbose treat
  if (verbose) cat( "\n" )

  ### net flows
  mfit[["gamma"]] <- colSums( sweep( mfit$pij , 1 , mfit$eta , "*" ) )
  mfit[["delta"]] <- N * ( mfit$gamma - mfit$eta )

  # return mcar
  if ( model == "MCAR" ) return( mfit )

  ### expected proportions under the model
  estimated.props <- this.expfun( c( mfit[["psi"]] , mfit[["rho"]] , mfit[["tau"]] ,mfit[["eta"]] , t( mfit[["pij"]] ) ) , k )

  # test value
  stopifnot( all.equal( sum( estimated.props ) , 1 ) )

  ### pseudo log-likelihood
  pll <- CountMatrix * log( estimated.props )
  mfit[["ll"]] <- sum( pll )

  ### unadjusted chi-distances

  # store observed counts
  mfit[["observed.counts"]] <- CountMatrix

  # store estimated counts
  mfit[["estimated.counts"]] <- estimated.props * N

  # calculate observed proportions
  observed.props <- CountMatrix / N

  # store observed counts
  mfit[["observed.props"]] <- observed.props

  # store estimated counts
  mfit[["estimated.props"]] <- estimated.props

  # store model info
  mfit[["model.info"]] <- c( n.counts , n.restr , n.parms )

  # store maxdiff
  mfit[["maxdiff"]] <- maxdiff

  # store fitting information
  if ( keep.info ) {
    mfit[["info.mat"]] <- info.mat[ !is.na( info.mat[,1] ) , ]
  }

  # return fit
  return( mfit )

}


