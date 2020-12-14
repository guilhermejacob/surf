# Wvec
modelC.WVec <- function( theta , CountMatrix ) {

  ##### object formatting

  # nummber of categories
  K <- dim( CountMatrix )[1] - 1

  # Obtain sample estimates of population flows as described in table 3.1 of Rojas et al. (2014)
  Nij <- as.matrix( CountMatrix[ -nrow( CountMatrix ) , -ncol( CountMatrix ) ] )
  Ri <- CountMatrix[ , ncol(CountMatrix) ][ - nrow( CountMatrix ) ]
  Cj <- CountMatrix[ nrow(CountMatrix) , ][ - ncol( CountMatrix ) ]
  M <- CountMatrix[ nrow( CountMatrix ) , ncol( CountMatrix ) ]
  N  <- sum( Nij ) + sum( Ri ) + sum( Cj ) + M

  # rebuild parameters form vector theta
  psi <- theta[ 1 ]
  rho <- theta[ 1+seq_len(K) ]
  tau <- theta[ (K+1) + seq_len(K) ]
  eta <- theta[ (2*K+1) + seq_len(K) ]
  pij <- theta[ 3*K+1 + seq_len(K^2) ]

  # test formats
  stopifnot( (3*K+1 + K^2 ) == length( theta) )

  # rebuild matrices
  pij <- matrix( pij , ncol = K , nrow = K , byrow = TRUE )

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )
  taunipij <- sweep( nipij , 1 , tau , "*" )
  taucnipij <- sweep( nipij , 1 , 1 - tau , "*" )
  tauni <- eta*tau
  taucni <- eta*(1-tau)
  taupij <- sweep( pij , 1 , tau , "*" )
  taucpij <- sweep( pij , 1 , 1 - tau , "*" )

  ##### estimating equations (Binder's W vector)

  # psi
  Wpsi <- ( sum( Nij ) + sum( Ri ) ) / psi - ( sum( Cj ) + M ) / ( 1- psi )

  # rho
  Wrho <- (rowSums( Nij ) / rho) - (Ri / ( 1 - rho ))

  # tau
  Wtau <- - rowSums( sweep( nipij , 2 , Cj / colSums( taucnipij ) , "*" ) ) + M * eta / sum( tauni )

  # eta
  Weta <-
    ( rowSums( Nij ) + Ri ) / eta +
    rowSums( sweep( taucpij , 2 , Cj / colSums( taucnipij ) , "*" ) ) +
    M * tau / sum( tauni ) - N

  # pij (unrestricted)
  Wpij <- Nij / pij
  Wpij <- sweep( Wpij , 1 , Ri / rowSums( pij ) , "+" )
  # Wpij <- sweep( Wpij , 1 , rowSums( outer( taucni , Cj / colSums( taucnipij ) ) ) , "+" )
  Wpij <- Wpij + outer( taucni , Cj / colSums( taucnipij ) )
  Wpij <- sweep( Wpij , 1 , M * tauni / sum( tauni ) , "+" )

  # lambda2 restriction
  lambda2 <-
    -( rowSums( Nij ) + Ri +
         rowSums( sweep( taucnipij , 2 , Cj / colSums( taucnipij ) , "*" ) ) +
         M * tauni / sum( tauni ) ) / N
  Wpij <- sweep( Wpij , 1 , N*lambda2 , "+" )

  # build Wvec
  c( Wpsi , Wrho , Wtau , Weta , t( Wpij ) )

}


# function for model variance calculation
modelC.variance <- function( xx , ww , res , design ) {

  # load objects
  Amat <- res[["observed.counts"]]
  K <- sqrt( prod( dim( Amat ) ) ) - 1
  this.theta <- c( unlist( res[ c( "psi" , "rho" , "tau" , "eta" ) ] ) , t( res[[ "pij" ]] ) )
  Nij <- res[["Nij"]]
  Kmat <- matrix( seq_len( prod( dim( Nij ) ) ) , nrow = nrow( Nij ) , byrow = TRUE )
  Ri <- res[["Ri"]]
  Cj <- res[["Cj"]]
  M <- res[["M"]]
  psi <- res[["psi"]]
  rho <- res[["rho"]]
  tau <- res[["tau"]]
  eta <- res[["eta"]]
  pij <- res[["pij"]]
  muij <- res[["muij"]]

  # yy array - see Rojas et al. (2014, p.294)
  yy <- array( 0  , dim = c( nrow( xx ) , nrow( Nij ) , ncol( xx ) ) )
  for ( r in seq_len( ncol( xx ) ) ) {
    kk <- stats::model.matrix( ~-1+. , data = xx[,r,drop = FALSE] , contrasts.arg = lapply( xx[,r, drop = FALSE ] , stats::contrasts, contrasts=FALSE ) , na.action = stats::na.pass )
    yy[ which( !is.na( xx[ , r ] ) ) , , r ] <- kk ; rm( kk )
  }

  # Create matrix of z variables - see Rojas et al. (2014, p.294)
  zz <- apply( xx , 2 , function(z) as.numeric( !is.na(z) ) )

  # Special variables - see Rojas et al. (2014, p.295)
  vv <- rowSums( yy[,,1] ) * rowSums( yy[,,2] ) + rowSums( yy[,,2] * (1 - zz[,1]) ) + rowSums( yy[,,1] * (1 - zz[,2]) ) + ( 1- zz[,1] ) * ( 1- zz[,2] )
  y1y2 <- array( 0  , dim = c( nrow( xx ) , nrow( Nij ) , ncol( Nij ) ) )
  for ( i in seq_len( nrow( Nij ) ) ) for ( j in seq_len( ncol( Nij ) ) ) y1y2[,i,j] <- yy[,i,1] * yy[,j,2]

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )
  taunipij <- sweep( nipij , 1 , tau , "*" )
  taucnipij <- sweep( nipij , 1 , 1 - tau , "*" )
  tauni <- eta*tau
  taucni <- eta*(1-tau)
  taupij <- sweep( pij , 1 , tau , "*" )
  taucpij <- sweep( pij , 1 , 1 - tau , "*" )

  ### psi

  # Calculate scores for estimating the variance of psi parameters
  u.psi <-
    ( rowSums( apply( y1y2 , c(1,3) , sum ) ) + rowSums( yy[,,1] ) * (1 - zz[,2] ) ) / psi -
    ( rowSums( yy[,,2] ) * (1 - zz[,1] ) + ( 1- zz[,1] ) * ( 1 - zz[,2] ) ) / ( 1 - psi )

  ### rho

  # Calculate scores for estimating the variance of rho parameters
  u.rho <- array( 0 , dim = c( nrow(xx) , nrow( Nij ) ) )
  for ( i in seq_len( nrow( Nij ) ) ) {
    u.rho[,i] <- yy[,i,1] * rowSums( yy[,,2] ) / rho[i] - yy[,i,1] * ( 1 - zz[,2] ) / ( 1 - rho[i] )
  }

  ### tau

  # Calculate scores for estimating the variance of tau parameters
  u.tau <- array( 0 , dim = c( nrow(xx) , nrow( Nij ) ) )
  for ( i in seq_len( nrow( Nij ) ) ) {
    u.tau[,i] <-
      ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * eta[i] / sum( tauni ) -
      rowSums( sweep( yy[,,2] * ( 1 - zz[,1] ) , 2 , nipij[i,] / colSums( taucnipij ) , "*" ) )
  }

  ### eta

  # Calculate scores for estimating the variance of eta parameters
  u.eta <- array( 0 , dim = c( nrow(xx) , nrow( Nij ) ) )
  for ( i in seq_len( nrow( Nij ) ) ) {
    u.eta[,i] <-
      ( ( yy[,i,1] * rowSums( yy[,,2] ) ) + ( yy[,i,1] * ( 1 - zz[,2] ) ) ) / eta[i] +
      rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , taucpij[i,] / colSums( taucnipij ) , "*" ) ) +
      ( 1- zz[,1] ) * ( 1 - zz[,2] ) * ( tau[i] / sum( tauni ) ) - 1
  }

  ### pij

  # Calculate scores for estimating the variance of pij parameters
  a.pij <- array( 0 , dim = c( nrow( xx ) , nrow( Nij ) , ncol( Nij ) ) )
  for ( i in seq_len( nrow( Nij ) ) ) for ( j in seq_len( ncol( Nij ) ) ) {
    a.pij[,i,j] <- ( y1y2[,i,j] / pij[i,j] ) +
      ( yy[,i,1] * ( 1 - zz[,2] ) ) +
      ( ( yy[,j,2] * ( 1 - zz[,1] ) ) * ( taucni[i] ) / colSums( taucnipij )[j] ) +
      ( ( ( 1 - zz[,1] ) * ( 1 - zz[,2] ) ) * ( ( tauni[i] ) / sum( tauni ) ) )
  }

  # lambda2 restriction
  lambda2 <-
    -( rowSums( Nij ) + Ri +
         rowSums( sweep( taucnipij , 2 , Cj / colSums( taucnipij ) , "*" ) ) +
         M * tauni / sum( tauni ) ) / N
  a.pij <- sweep( a.pij , 2 , lambda2 , "+" )

  # coerce to matrix
  u.pij <- matrix( 0 , nrow = dim( a.pij )[1] , ncol = K^2 , byrow = TRUE )
  for ( i in seq_len( nrow( Nij ) ) ) u.pij[,Kmat[i,]] <- a.pij[,i,]

  ### matrix of linearized variables

  # build Umat
  Umat <- do.call( cbind , list( u.psi , u.rho , u.tau , u.eta , u.pij ) )

  # test equality (within some tolerance)
  stopifnot( all.equal( colSums( Umat * ww ) , modelC.WVec( this.theta , Amat ) , check.attributes = FALSE ) )

  ### calculate jacobian

  # jacobian matrix
  Jmat <- numDeriv::jacobian( modelC.WVec , this.theta , method = "complex" , side = NULL , CountMatrix = Amat )

  # inverse of the jacobian matrix
  Jmat.inv <- MASS::ginv( Jmat )

  # # calculate variance #1
  # vmat0 <- survey::svyrecvar( sweep( Umat , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  # full.vmat0 <- Jmat.inv %*% vmat0 %*% t( Jmat.inv )

  # calculate variance #2
  Umat <- t( apply( Umat , 1 , function(z) crossprod( -t( Jmat.inv ) , z ) ) )
  u.psi <- Umat[ , 1 ]
  u.rho <- Umat[ , 2 ]
  u.tau <- Umat[ , 3 ]
  u.eta <- Umat[ , seq( 4 , 4 + K - 1 ) ]
  u.pij <- Umat[ , seq( 4 + K , 4 + K + ( K^2 - 1 ) ) ]
  a.pij <- array( 0 , dim = c( nrow( xx ) , nrow( Nij ) , ncol( Nij ) ) )
  for ( j in seq_len( ncol( Nij ) ) ) {
    a.pij[,,j] <- u.pij[ , Kmat[ j , ] ]
  }
  full.vmat <- survey::svyrecvar( sweep( Umat , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

  ##### other variances

  # net flows
  u.nipij <- array( 0 , dim = c( nrow( xx ) , nrow( Nij ) , ncol( Nij ) ) )
  for ( i in seq_len( nrow(Nij) ) ) for ( j in seq_len( ncol( Nij ) ) ) {
    u.nipij[,i,j] <- ( pij[i,j] * u.eta[,i] + eta[i] * a.pij[,i,j] )
  }
  # u.nipij <- matrix( u.nipij , nrow = dim( u.pij )[1] , byrow = FALSE )
  # nipij.vmat <- survey::svyrecvar( sweep( u.nipij , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  # nipij.vmat <- matrix( diag( nipij.vmat ) , nrow = nrow( Nij ) , byrow = FALSE )

  # gross flows
  a.muij <- array( 0 , dim = c( nrow( xx ) , nrow( Nij ) , ncol( Nij ) ) )
  for ( i in seq_len( nrow(Nij) ) ) for ( j in seq_len( ncol( Nij ) ) ) {
    a.muij[,i,j] <- nipij[i,j] + N * u.nipij[,i,j]
  }
  u.muij <- matrix( 0 , nrow = dim( a.muij )[1] , ncol = K^2 , byrow = TRUE )
  for ( i in seq_len( nrow( Nij ) ) ) u.muij[,Kmat[i,]] <- a.muij[,i,]
  muij.vmat <- survey::svyrecvar( sweep( u.muij , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )
  rm( a.muij )

  # final distribution
  u.gamma <- apply( u.nipij , c(1,3) , sum )
  gamma.vmat <- survey::svyrecvar( sweep( u.gamma , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

  # delta
  delta <- N * ( colSums( nipij ) - eta )
  u.delta <- sweep( N * ( u.gamma - u.eta ) , 2 , res[["gamma"]] - res[["eta"]] , "+" )
  delta.vmat <- survey::svyrecvar( sweep( u.delta , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata )

  ##### split full matrix

  # non-response
  psi.vmat <- diag( full.vmat )[1]
  rho.vmat <- diag( full.vmat )[2]
  tau.vmat <- diag( full.vmat )[3]

  # unobserved process
  eta.vmat <- full.vmat[ seq( 4 , 4 + K - 1 ) , seq( 4 , 4 + K - 1 ) ]
  pij.vmat <- full.vmat[ seq( 4 + K , 4 + K + ( K^2 - 1 ) ) , seq( 4 + K , 4 + K + ( K^2 - 1 ) ) ]

  # collect variances from block diagonal matrix
  muij.vmat <- matrix( diag( muij.vmat ) , nrow = nrow( Nij ) , byrow = FALSE )
  pij.vmat <- matrix( diag( pij.vmat ) , nrow = nrow( Nij ) , byrow = FALSE )

  # build list of variances
  mvar <-
    list(
      "psi" = matrix( psi.vmat ) ,
      "rho" = matrix( rho.vmat ) ,
      "tau" = matrix( tau.vmat ) ,
      "eta" = eta.vmat ,
      "pij" = pij.vmat ,
      "muij" = muij.vmat ,
      "gamma" = gamma.vmat ,
      "delta" = delta.vmat )

  # return list of variances
  return( mvar )

}


# initial values
modelC.initial <- function( CountMatrix ) {

  # non-response cells
  Nij <- CountMatrix[ -nrow( CountMatrix ) , -ncol( CountMatrix ) ]
  Ri <- CountMatrix[ , ncol(CountMatrix) ][ - nrow( CountMatrix ) ]
  Cj <- CountMatrix[ nrow( CountMatrix ) , ][ - ncol( CountMatrix ) ]
  M <- CountMatrix[ nrow( CountMatrix ) , ncol( CountMatrix ) ]
  N <- sum( CountMatrix )

  # Obtain maximum pseudo-likelihood estimates for response model parameters
  # rho, and tau according to Result 4.12 of Rojas (2014, p.54-55)
  psi <- ( sum( Nij ) + sum( Ri ) ) / N
  rho <- rowSums( Nij ) / ( rowSums( Nij ) + Ri )

  # Obtain starting values for estimating superpopulation model flow parameters
  # psi, eta and pij according to Result 4.13 of Rojas (2014, p.62)
  tau <- rep( M / ( sum( Cj ) + M ) , nrow( Nij ) )
  eta <- rowSums( Nij ) / sum( Nij )
  pij <- sweep( Nij , 1 , rowSums( Nij ) , "/" )

  # combine initial values
  c( psi , rho , tau , eta , t( pij ) )

}

# model fitting
modelC.fitting <- function( theta , CountMatrix , tol = 1/sum( CountMatrix ) , maxit = 500 , verbose = FALSE ) {

  # prepare iterations
  maxdiff <- 1
  v <- 0
  ndig <- ceiling( -log(tol,10) ) + 2
  last.loglik <- -Inf
  # last.loglik <- modelC.loglik( theta , CountMatrix )

  # process tracker
  if (verbose) cat( "\n" )

  # iterations
  while ( v < maxit ) {
    v = v+1
    theta <- matrix( theta , ncol = 1 )
    Wvec <- matrix( modelC.WVec( theta , CountMatrix ) , ncol = 1 )
    Jmat <- numDeriv::jacobian( modelC.WVec , theta , method = "complex" , side = NULL , CountMatrix = CountMatrix )
    Jmat.inv <- MASS::ginv( Jmat , tol = 1e-10 )
    this.theta <- theta - crossprod( Jmat.inv , Wvec )
    # while ( any( this.theta > 1 | this.theta < 0 ) ) this.theta[ this.theta > 1 | this.theta < 0 ] <- (( theta + this.theta ) / 2) [ this.theta > 1 | this.theta < 0 ]
    delta.theta <- c( this.theta - theta )
    maxdiff<- max( abs( delta.theta ) )
    this.loglik <- modelC.loglik( this.theta , CountMatrix )
    delta.loglik <- this.loglik - last.loglik

    # process tracker
    if (verbose) cat( sprintf( paste0( "iteration: %5s\t maxdiff: %1." , ndig , "f \t (p)ll: %-10." , ndig , "f \n" ) , v , maxdiff , this.loglik ) )

    # break loop
    if ( maxdiff < tol ) break()

    # else change
    theta <- this.theta
    last.loglik <- this.loglik

  }

  # process tracker
  if (verbose) cat( "\n" )

  c( this.theta )

}

# pseudo-likelihood
modelC.loglik <- function( theta , CountMatrix ) {

  # nummber of categories
  K <- dim( CountMatrix )[1] - 1

  # rebuild parameters form vector theta
  psi <- theta[ 1 ]
  rho <- theta[ 1+seq_len(K) ]
  tau <- theta[ (K+1) + seq_len(K) ]
  eta <- theta[ (2*K+1) + seq_len(K) ]
  pij <- theta[ 3*K+1 + seq_len(K^2) ]

  # test formats
  stopifnot( (3*K+1 + K^2 ) == length( theta) )

  # rebuild matrices
  pij <- matrix( pij , ncol = K , nrow = K , byrow = TRUE )

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )
  rhonipij <- sweep( nipij , 1 , rho , "*" )
  rhocnipij <- sweep( nipij , 1 , 1 - rho , "*" )
  taunipij <- sweep( nipij , 1 , tau , "*" )
  taucnipij <- sweep( nipij , 1 , 1 - tau , "*" )

  # matrix blocks
  Part.Nij <- psi * rhonipij
  Part.Cj <- colSums( (1 - psi) * taucnipij )
  Part.Ri <- rowSums( psi * rhocnipij )
  Part.M <- sum( (1 - psi ) * taunipij )

  # build matrix
  expected.props <- rbind( cbind( Part.Nij , Part.Cj ) , c( Part.Ri , Part.M ) )

  # evaluate
  sum( CountMatrix * log( expected.props ) ) # unconstrained

}

# expected proportions
modelC.expected <- function( theta , K ) {

  # rebuild parameters form vector theta
  psi <- theta[ 1 ]
  rho <- theta[ 1+seq_len(K) ]
  tau <- theta[ (K+1) + seq_len(K) ]
  eta <- theta[ (2*K+1) + seq_len(K) ]
  pij <- theta[ 3*K+1 + seq_len(K^2) ]

  # test formats
  stopifnot( (3*K+1 + K^2 ) == length( theta) )

  # rebuild matrices
  pij <- matrix( pij , ncol = K , nrow = K , byrow = TRUE )

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )
  rhonipij <- sweep( nipij , 1 , rho , "*" )
  rhocnipij <- sweep( nipij , 1 , 1 - rho , "*" )
  taunipij <- sweep( nipij , 1 , tau , "*" )
  taucnipij <- sweep( nipij , 1 , 1 - tau , "*" )

  # matrix blocks
  Part.Nij <- psi * rhonipij
  Part.Cj <- colSums( (1 - psi) * taucnipij )
  Part.Ri <- rowSums( psi * rhocnipij )
  Part.M <- sum( (1 - psi ) * taunipij )

  # build matrix
  rbind( cbind( Part.Nij , Part.Cj ) , c( Part.Ri , Part.M ) )

}
