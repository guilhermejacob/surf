# function for model fitting
ipf <- function( CountMatrix , model , tol = NULL , maxit = 500 , verbose = FALSE ) {

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

  # initial settings for iterations
  maxdiff <- 1
  v <- 0
  this.step <- 1 # for tau in model D
  last.likelihood <- 0

  # define log-likelihood
  this.loglik <- switch( model , A = pll.modA , B = pll.modB , C = pll.modC , D = pll.modD )

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
      nipij <- sweep( pijv , 1 , etav , "*" )

      # calculate eta v+1
      etav <-
        ( rowSums( Nij ) + Ri + rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , Cj , "*" ) ) ) /
        ( sum( Nij ) + sum( Ri ) + sum( Cj ) )

      # calculate pij v+1
      pijv <-
        sweep( Nij + sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , Cj , "*" ) , 1 ,
               rowSums( Nij ) + rowSums( sweep( sweep( nipij , 2 , colSums( nipij ) , "/" ) , 2 , Cj , "*" ) ) , "/" )

      # calculate differences
      these.diffs <- c( etav - eta0 , pijv - pij0 )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these.diffs ) )

      # process tracker
      if (verbose) cat( sprintf( paste0( "iteration: %5s \t maxdiff: %1." , ceiling( -log(tol,10) ) + 2 , "f \r" ) , v , maxdiff ) )

      # test for non-convergence
      if ( v >= maxit ) {
        stop( "Algorithm has not converged after " , v , " iterations." )
        break()
      }

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
        "Nij" = Nij ,
        "Ri" = Ri ,
        "Cj" = Cj ,
        "M" = M ,
        "psi" = psi ,
        "rho" = rho ,
        "tau" = tau ,
        "eta" = etav ,
        "pij" = pijv ,
        "muij" = N * sweep( pijv , 1 , etav , "*" ) ,
        "gamma" = colSums( sweep( pijv , 1 , etav , "*" ) ) )
  } , B = {
    # Obtain maximum pseudo-likelihood estimates for response model parameters
    # rho, and tau according to Result 4.7 of Rojas (2014, p.46-47)
    rho <- sum( Nij ) / ( sum( Nij ) + sum( Ri ) )
    tau <- M / ( sum( Cj ) + M )

    # Obtain starting values for estimating superpopulation model flow parameters
    # psi, eta and pij according to Result 4.8 of Rojas (2014, p.49)
    # psi0 <- psiv <- rep( ( sum( Nij ) + sum( Ri ) ) / N , nrow( Nij ) ) # original
    psi0 <- psiv <- ( rowSums( Nij ) + Ri ) / N # modified #1
    # psi0 <- psiv <- ( rowSums( Nij ) + Ri ) / ( sum( Nij ) + sum( Ri ) ) # modified #2
    eta0 <- etav <- rowSums( Nij ) / sum( Nij )
    pij0 <- pijv <- sweep( Nij , 1 , rowSums( Nij ) , "/" )

    # Obtain maximum pseudo-likelihood estimates for superpopulation model flow parameters
    # psi, eta and pij according to Result 4.8 of Rojas (2014, p.296)
    while( maxdiff >= tol ) {

      # count iteration
      v <- v+1

      # calculate auxiliary stats
      nipij <- sweep( pijv , 1 , etav , "*" )
      psicni <- ( 1 - psiv ) * etav
      psicnipij <- sweep( nipij , 1 , 1 - psiv , "*" )

      # calculate psi
      psiv <-
        ( rowSums( Nij ) + Ri ) /
        ( rowSums( Nij ) + Ri +
            rowSums( sweep( psicnipij , 2 , Cj / colSums( psicnipij ) , "*" ) ) +
            ( M * ( psicni / sum( psicni ) ) ) )

      # calculate eta v+1
      etav <-
        ( rowSums( Nij ) + Ri +
            rowSums( sweep( psicnipij , 2 , Cj / colSums( psicnipij ) , "*" ) ) +
            ( M * ( psicni / sum( psicni ) ) ) ) / N

      # calculate pij v+1
      pijv <-
        sweep( Nij + sweep( psicnipij , 2 , Cj / colSums( psicnipij ) , "*" ) , 1 ,
               rowSums( Nij + sweep( psicnipij , 2 , Cj / colSums( psicnipij ) , "*" ) ) , "/" )

      # # check value consistency #3
      # psiv <- ( psiv + psi0 ) / 2
      # etav <- ( etav + eta0 ) / 2
      # pijv <- ( pijv + pij0 ) / 2
      # while ( max( abs( psiv - psi0 ) ) > this.step ) psiv <- ( psiv + psi0 ) / 2
      # this.step <- max( abs( psiv - psi0 ) )

      # calculate differences
      these.diffs <- c( c( psiv - psi0 ) , c( etav - eta0 ) , c( pijv - pij0 ) )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these.diffs ) )

      # process tracker
      if (verbose) cat( sprintf( paste0( "iteration: %5s \t maxdiff: %1." , ceiling( -log(tol,10) ) + 2 , "f \r" ) , v , maxdiff ) )

      # test for non-convergence
      if ( v >= maxit ) {
        stop( "Algorithm has not converged after " , v , " iterations." )
        break()
      }

      # store estimates for next iteration
      psi0 <- psiv
      eta0 <- etav
      pij0 <- pijv

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
        "muij" = N * sweep( pijv , 1 , etav , "*" ) ,
        "gamma" = colSums( sweep( pijv , 1 , etav , "*" ) ) )

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

      # process tracker
      if (verbose) cat( sprintf( paste0( "iteration: %5s \t maxdiff: %1." , ceiling( -log(tol,10) ) + 2 , "f \r" ) , v , maxdiff ) )

      # test for non-convergence
      if ( v >= maxit ) {
        stop( "Algorithm has not converged after " , v , " iterations." )
        break()
      }

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
        "Nij" = Nij ,
        "Ri" = Ri ,
        "Cj" = Cj ,
        "M" = M ,
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
    psi <- ( sum( Nij ) + sum( Ri ) ) / ( sum( Nij ) + sum( Ri ) + sum( Cj ) +M )

    # Obtain starting values for estimating superpopulation model flow parameters
    # rhoj, tauj, eta and pij
    rhoj0 <- rhojv <- rep( sum( Nij ) / ( sum( Nij ) + sum( Ri ) ) , ncol( Nij ) )
    # tauj0 <- taujv <- rep( M / ( sum( Cj ) + M ) , ncol( Nij ) )
    tauj0 <- taujv <- ( M / ( Cj + M ) )
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

      # calculate rhoj
      rhojv <-
        colSums( Nij ) /
        ( colSums( Nij ) + colSums( sweep( sweep( nipijrhoc , 1 , Ri , "*" ) , 1 , rowSums( nipijrhoc ) , "/" ) ) )

      # calculate tauj #2
      taujv <- 1 - ( Cj * sum( nipijtau ) ) / ( M * colSums( nipij ) ) # Stasny
      # taujv <- ( M * colSums( nipij ) - Cj * ( sum( nipijtau ) - colSums( nipijtau ) ) ) / ( ( M + Cj ) * colSums( nipij ) ) # Hanwen & Andrés

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

      # # check value consistency #1
      # while ( any( c( rhojv , taujv ) < 0 | c( rhojv , taujv ) > 1 ) ) {
      #   rhojv <- ( rhojv + rhoj0 ) / 2
      #   taujv <- ( taujv + tauj0 ) / 2
      #   etav <- ( etav + eta0 ) / 2
      #   pijv <- ( pijv + pij0 ) / 2
      # }

      # check value consistency #2
      while ( any( rhojv < 0 | rhojv > 1 ) ) rhojv <- ( rhojv + rhoj0 ) / 2
      while ( any( taujv < 0 | taujv > 1 ) ) taujv <- ( taujv + tauj0 ) / 2

      # # check value consistency #3
      # while ( max( abs( taujv - tauj0 ) ) > this.step ) taujv <- ( taujv + tauj0 ) / 2
      # this.step <- max( abs( taujv - tauj0 ) )

      # damping
      rhojv <- ( rhojv + rhoj0 ) / 2
      taujv <- ( taujv + tauj0 ) / 2
      etav <- ( etav + eta0 ) / 2
      pijv <- ( pijv + pij0 ) / 2

      # calculate differences
      these.diffs <- c( c( rhojv - rhoj0 ) , c( taujv - tauj0 ) , c( etav - eta0 ) , c( pijv - pij0 ) )

      # calculate maximum absolute difference
      maxdiff <- max( abs( these.diffs ) )

      # calculate pseudo-likelihood
      iter.likelihood <- this.loglik( psi , rhojv , taujv , etav , pijv , CountMatrix )
      delta.likelihood <- iter.likelihood - last.likelihood

      # process tracker
      if (verbose) cat( sprintf( paste0( "iteration: %5s \t maxdiff: %1." , ceiling( -log(tol,10) ) + 2 , "f \r" ) , v , maxdiff ) )

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
      last.likelihood <- iter.likelihood

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
        "muij" = N * sweep( pijv , 1 , etav , "*" ) ,
        "gamma" = colSums( sweep( pijv , 1 , etav , "*" ) ) )

  } )

  # verbose treat
  if (verbose) cat( "\n" )

  ### net flows
  mfit[["delta"]] <- N * ( mfit$gamma - mfit$eta )

  # return mcar
  if ( model == "MCAR" ) return( mfit )

  ### expected proportions under the model
  estimated.props <- if ( model == "A" ) {
    cbind( rbind( sweep( res[["pij"]] , 1 , res[["eta"]] * res[["psi"]] , "*" ) * res[["rho"]] ,
                  colSums( sweep( res[["pij"]] , 1 , res[["eta"]] * ( 1 - res[["psi"]] ) , "*" ) * ( 1 - res[["tau"]] ) ) ) ,
           c( rowSums( sweep( res[["pij"]] , 1 , res[["eta"]] * res[["psi"]] , "*" ) * ( 1 - res[["rho"]] ) ) ,
              sum( sweep( res[["pij"]] , 1 , res[["eta"]] * ( 1 - res[["psi"]] ) , "*" ) *  res[["tau"]] ) ) )
  } else if ( model == "B" ) {
    cbind( rbind( sweep( res[["pij"]] , 1 , res[["eta"]] * res[["psi"]] , "*" ) * res[["rho"]] ,
                  colSums( sweep( res[["pij"]] , 1 , res[["eta"]] * ( 1 - res[["psi"]] ) , "*" ) * ( 1 - res[["tau"]] ) ) ) ,
           c( rowSums( sweep( res[["pij"]] , 1 , res[["eta"]] * res[["psi"]] , "*" ) * ( 1 - res[["rho"]] ) ) ,
              sum( sweep( res[["pij"]] , 1 , res[["eta"]] * ( 1 - res[["psi"]] ) , "*" ) *  res[["tau"]] ) ) )
  } else if ( model == "C" ) {
    cbind( rbind( sweep( sweep( res[["pij"]] , 1 , res[["eta"]] * res[["psi"]] , "*" ) , 1 , res[["rho"]] , "*" ) ,
                  colSums( sweep( sweep( res[["pij"]] , 1 , res[["eta"]] * ( 1 - res[["psi"]] ) , "*" ) , 1 , 1 - res[["tau"]] , "*" ) ) ) ,
           c( rowSums( sweep( sweep( res[["pij"]] , 1 , res[["eta"]] * res[["psi"]] , "*" ) , 1 , 1 - res[["rho"]] , "*" ) ) ,
              sum( sweep( sweep( res[["pij"]] , 1 , res[["eta"]] * ( 1 - res[["psi"]] ) , "*" ) , 1 , res[["tau"]] , "*" ) ) ) )
  } else if ( model == "D" ) {
    cbind( rbind( sweep( sweep( res[["pij"]] , 1 , res[["eta"]] * res[["psi"]] , "*" ) , 2 , res[["rho"]] , "*" ) ,
                  colSums( sweep( sweep( res[["pij"]] , 1 , res[["eta"]] * ( 1 - res[["psi"]] ) , "*" ) , 2 , 1 - res[["tau"]] , "*" ) ) ) ,
           c( rowSums( sweep( sweep( res[["pij"]] , 1 , res[["eta"]] * res[["psi"]] , "*" ) , 2 , 1 - res[["rho"]] , "*" ) ) ,
              sum( sweep( sweep( res[["pij"]] , 1 , res[["eta"]] * ( 1 - res[["psi"]] ) , "*" ) , 2 , res[["tau"]] , "*" ) ) ) )
  }
  dimnames( estimated.props ) <- dimnames( CountMatrix )

  # test value
  stopifnot( all.equal( sum( estimated.props ) , 1 ) )

  ### log pseudo-likelihood
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

  # chi-distance matrix
  chimat <- ( observed.props - estimated.props )^2 / estimated.props
  mfit[["chimat"]] <- chimat

  # calculate unadjusted test score
  chiscore <- N * sum( chimat )

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

  # store model info
  mfit[["model.info"]] <- c( n.counts , n.restr , n.parms )

  # store maxdiff
  mfit[["maxdiff"]] <- maxdiff

  # return fit
  return( mfit )

}

# pseudo-likelihood function
pll.modA <- function( psi , rho , tau , eta , pij , CountVector ) {

  # rebuild matrices
  pij <- matrix( pij , nrow = sqrt( length( pij ) ) , byrow = TRUE )
  CountMatrix <- matrix( CountVector , nrow = sqrt( length( CountVector ) ) , byrow = TRUE )

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )

  # matrix blocks
  Part.Nij <- sweep( nipij , 1 , psi , "*" ) * rho
  Part.Cj <- colSums( nipij * ( 1 - psi ) * ( 1 - tau ) )
  Part.Ri <- rowSums( nipij * psi * ( 1 - rho ) )
  Part.M <- sum( nipij * ( 1 - psi ) *  tau )

  # build matrix
  expected.props <- rbind( cbind( Part.Nij , Part.Cj ) , c( Part.Ri , Part.M ) )

  # evaluate
  sum( CountMatrix * expected.props )

}

