variance_fun <- function( xx , ww , res , design ) {

  switch( res$model ,
          A = { variance.fun <- modelA.variance } ,
          B = { variance.fun <- modelB.variance } ,
          C = { variance.fun <- modelC.variance } ,
          D = { variance.fun <- modelD.variance } )

  variance.fun( xx , ww , res , design )

}

fitting_fun <- function( CountMatrix , model , tol = max( 1e-6 , 1/sum( CountMatrix ) ) , maxit = 500 , verbose = FALSE , keep.info = FALSE ) {

  # define function to estimate initial values for the fitting procedure
  switch( model ,
          A = { startval.fun <- modelA.initial } ,
          B = { startval.fun <- modelB.initial } ,
          C = { startval.fun <- modelC.initial } ,
          D = { startval.fun <- modelD.initial } )

  # define (pseudo) log-likelihood function
  switch( model ,
          A = { loglik.fun <- modelA.loglik } ,
          B = { loglik.fun <- modelB.loglik } ,
          C = { loglik.fun <- modelC.loglik } ,
          D = { loglik.fun <- modelD.loglik } )

  # define function to estimate Binder's W vector
  switch( model ,
          A = { Wvec.fun <- modelA.WVec } ,
          B = { Wvec.fun <- modelB.WVec } ,
          C = { Wvec.fun <- modelC.WVec } ,
          D = { Wvec.fun <- modelD.WVec } )

  # define function to estimate Binder's W vector
  switch( model ,
          A = { expect.fun <- modelA.expected } ,
          B = { expect.fun <- modelB.expected } ,
          C = { expect.fun <- modelC.expected } ,
          D = { expect.fun <- modelD.expected } )

  # number of categories
  K <- sqrt( prod( dim( CountMatrix ) ) ) - 1

  # number of model parameters
  switch( model ,
          A = { n.parms <- K^2 +   K + 3 } ,
          B = { n.parms <- K^2 + 2*K + 2 } ,
          C = { n.parms <- K^2 + 3*K + 1 } ,
          D = { n.parms <- K^2 + 3*K + 1 } )

  ### generalized fitting

  # create algorithm info
  if ( keep.info ) {
    info.mat <- matrix( as.numeric( NA ) , ncol = n.parms + 3 , nrow = maxit )
  }

  # prepare iterations
  maxdiff <- Inf
  v <- 0
  ndig <- ceiling( -log(tol,10) ) + 2
  last.loglik <- -Inf

  # calculate initial values
  theta <- startval.fun( CountMatrix )

  # iterations
  while ( maxdiff > tol ) {

    # calculate new parameters
    v = v+1
    theta <- matrix( theta , ncol = 1 )
    Wvec <- matrix( Wvec.fun( theta , CountMatrix ) , ncol = 1 )
    Jmat <- numDeriv::jacobian( Wvec.fun , theta , method = "complex" , side = NULL , CountMatrix = CountMatrix )
    Jmat.inv <- MASS::ginv( Jmat )
    this.theta <- theta - crossprod( Jmat.inv , Wvec )

    # damping
    # while ( any( this.theta > 1 | this.theta < 0 ) ) this.theta[ this.theta > 1 | this.theta < 0 ] <- (( theta + this.theta ) / 2) [ this.theta > 1 | this.theta < 0 ]
    while ( any( this.theta > 1 | this.theta < 0 ) ) this.theta <- (( theta + this.theta ) / 2)

    # results
    delta.theta <- c( this.theta - theta )
    maxdiff<- max( abs( delta.theta ) )
    this.loglik <- loglik.fun( this.theta , CountMatrix )
    delta.loglik <- this.loglik - last.loglik

    # process tracker
    if (verbose) cat( sprintf( paste0( "iteration: %5s\t maxdiff: %1." , ndig , "f \t (p)ll: %-10." , ndig , "f \n" ) , v , maxdiff , this.loglik ) )

    # store fitting information
    if ( keep.info ) {
      info.mat[ v , ] <- c( v , theta , maxdiff , this.loglik )
    }

    # test for non-convergence
    if ( v >= maxit ) {
      stop( "Algorithm has not converged after " , v , " iterations." )
      break()
    }

    # break loop
    if ( maxdiff < tol ) break()

    # else change
    theta <- this.theta
    last.loglik <- this.loglik

  }

  # process tracker
  if (verbose) cat( "\n" )

  return( this.theta )

  ### format results

  # Obtain sample estimates of population flows as described in table 3.1 of Rojas et al. (2014)
  Nij <- as.matrix( CountMatrix[ -nrow( CountMatrix ) , -ncol( CountMatrix ) ] )
  Ri <- CountMatrix[ , ncol(CountMatrix) ][ - nrow( CountMatrix ) ]
  Cj <- CountMatrix[ nrow(CountMatrix) , ][ - ncol( CountMatrix ) ]
  M <- CountMatrix[ nrow( CountMatrix ) , ncol( CountMatrix ) ]
  N  <- sum( Nij ) + sum( Ri ) + sum( Cj ) + M

  # rebuild parameters form vector theta
  switch( model , A = {

    psi <- theta[ 1 ]
    rho <- theta[ 2 ]
    tau <- theta[ 3 ]
    eta <- theta[ 3 + seq_len( K ) ]
    pij <- theta[ seq( 4 + K , ( 4 + K ) + K^2 - 1 ) ]

  } , B = {

    psi <- theta[ 1:K ]
    rho <- theta[ K+1 ]
    tau <- theta[ K+2 ]
    eta <- theta[ seq( K+3 , K+3 + ( K - 1 ) ) ]
    pij <- theta[ seq( 2*K+3 , ( 2*K+3 ) + ( K^2 - 1 ) ) ]

  } , C = {
    psi <- this.theta[ 1 ]
    rho <- this.theta[ 1+seq_len(K) ]
    tau <- this.theta[ (K+1) + seq_len(K) ]
    eta <- this.theta[ (2*K+1) + seq_len(K) ]
    pij <- this.theta[ 3*K+1 + seq_len(K^2) ]
  } , D = {
    psi <- this.theta[ 1 ]
    rho <- this.theta[ 1+seq_len(K) ]
    tau <- this.theta[ (K+1) + seq_len(K) ]
    eta <- this.theta[ (2*K+1) + seq_len(K) ]
    pij <- this.theta[ 3*K+1 + seq_len(K^2) ]
  } )

  # rebuild matrices
  pij <- matrix( pij , ncol = K , nrow = K , byrow = TRUE )

  # other estimates
  nipij <- sweep( pij , 1 , eta , "*" )
  muij <- N * nipij
  gamma <- colSums( nipij )
  delta <- N * ( gamma - eta )

  # return final estimates
  mfit <-
    list( "model" = model ,
          "iter" = v ,
          "N" = N ,
          "Nij" = Nij ,
          "Ri" = Ri ,
          "Cj" = Cj ,
          "M" = M ,
          "psi" = psi ,
          "rho" = rho ,
          "tau" = tau ,
          "eta" = eta ,
          "pij" = pij ,
          "muij" = muij ,
          "gamma" = gamma ,
          "delta" = delta )

  # add (pseudo) log-likelihood
  mfit[["ll"]] <- this.loglik

  # add expected proportions and counts
  expected.props <- expect.fun( theta , K = K )

  # store observed counts
  mfit[["observed.counts"]] <- CountMatrix
  mfit[["observed.props"]] <- CountMatrix / N

  # store estimated counts
  mfit[["estimated.counts"]] <- N * expected.props
  mfit[["estimated.props"]] <- expected.props

  # store maxdiff
  mfit[["maxdiff"]] <- maxdiff

  # return fitted model information
  return( mfit )

}
