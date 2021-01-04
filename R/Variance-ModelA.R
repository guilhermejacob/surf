# Wvec
modelA.WVec <- function( theta , CountMatrix ) {

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
  rho <- theta[ 2 ]
  tau <- theta[ 3 ]
  eta <- theta[ 3 + seq_len( K ) ]
  pij <- theta[ (K+3) + seq_len(K^2) ]

  # rebuild matrices
  pij <- matrix( pij , ncol = K , nrow = K , byrow = TRUE )

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )

  ##### estimating equations (Binder's W vector)

  # psi
  Wpsi <- (( sum( Nij ) + sum( Ri ) ) / psi) - (( sum( Cj ) + M ) / ( 1- psi ))


  # rho
  Wrho <- (sum( Nij ) / rho) - (sum( Ri ) / ( 1 - rho ))

  # tau
  Wtau <- - (sum( Cj ) / ( 1 - tau )) + (M / tau)

  # eta
  Weta <-
    ( rowSums( Nij ) + Ri ) / eta +
    rowSums( sweep( pij , 2 , Cj / colSums( nipij ) , "*" ) ) +
    M - N

  # pij (unrestricted)
  Wpij <- (Nij / pij)
  Wpij <- sweep( Wpij , 1 , Ri / rowSums( pij ) , "+" )
  Wpij <- Wpij + outer( eta , Cj / colSums( nipij ) )
  Wpij <- sweep( Wpij , 1 , M * ( eta / sum( nipij ) ) , "+" )

  # pij = zero
  Wpij[ is.na( Wpij ) ] <- 0

  # lambda2 restriction
  lambda2 <-
    -( rowSums( Nij ) + Ri +
         rowSums( sweep( nipij , 2 , Cj / colSums( nipij ) , "*" ) ) +
         M * eta ) / N
  Wpij <- sweep( Wpij , 1 , N*lambda2 , "+" )

  # build Wvec
  c( Wpsi , Wrho , Wtau , Weta , t( Wpij ) )

}

# function for model variance calculation
modelA.linearization <- function( xx , ww , res , design ) {

  # load objects
  Amat <- res[["observed.counts"]]
  K <- sqrt( prod( dim( Amat ) ) ) - 1
  this.theta <- c( unlist( res[ c( "psi" , "rho" , "tau" , "eta" ) ] ) , t( res[[ "pij" ]] ) )
  Nij <- res[["Nij"]]
  Ri <- res[["Ri"]]
  Cj <- res[["Cj"]]
  M <- res[["M"]]
  N <- sum( Amat )
  psi <- res[["psi"]]
  rho <- res[["rho"]]
  tau <- res[["tau"]]
  eta <- res[["eta"]]
  pij <- res[["pij"]]
  muij <- res[["muij"]]
  Kmat <- matrix( seq_len( prod( dim( Nij ) ) ) , nrow = nrow( Nij ) , byrow = TRUE )

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

  # calculate auxiliary stats
  nipij <- sweep( pij , 1 , eta , "*" )

  ### psi

  # Calculate scores for estimating the variance of psi parameters
  u.psi <-
    (( apply( y1y2 , 1 , sum ) + rowSums( yy[,,1] ) * (1 - zz[,2] ) ) / psi) -
    (( rowSums( yy[,,2] ) * (1 - zz[,1] ) + ( 1- zz[,1] ) * ( 1 - zz[,2] ) ) / ( 1 - psi ))

  ### rho

  # Calculate scores for estimating the variance of rho parameters
  u.rho <- ( apply( y1y2 , 1 , sum ) / rho ) - ( ( rowSums( yy[,,1] ) * ( 1 - zz[,2]) ) / ( 1 - rho ) )

  ### tau

  # Calculate scores for estimating the variance of tau parameters
  u.tau <- ( ( 1 - zz[,1] ) * ( 1 - zz[,2] ) / tau ) - ( ( rowSums( yy[,,2] ) * ( 1 - zz[,1] ) ) / ( 1 - tau ) )

  ### eta

  # Calculate scores for estimating the variance of eta parameters
  u.eta <- array( 0 , dim = c( nrow(xx) , nrow( Nij ) ) )
  for ( i in seq_len( nrow( Nij ) ) ) {
    u.eta[,i] <-
      rowSums( y1y2[,i,] ) / eta[i] +
      ( yy[,i,1] * ( 1 - zz[,2] ) ) / eta[i] +
      # rowSums( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , pij[i,] / colSums( nipij ) , "*" ) ) +
      rowSums( sweep( sweep( yy[,,2] * (1 - zz[,1] ) , 2 , pij[i,] , "*" ) , 2 , colSums( nipij ) , "/" ) ) +
      ( 1- zz[,1] ) * ( 1 - zz[,2] ) - 1
  }

  ### pij

  # Calculate scores for estimating the variance of pij parameters
  a.pij <- array( 0 , dim = c( nrow( xx ) , nrow( Nij ) , ncol( Nij ) ) )
  for ( i in seq_len( nrow( Nij ) ) ) for ( j in seq_len( ncol( Nij ) ) ) {
    a.pij[,i,j] <-
      ( y1y2[,i,j] / pij[i,j] ) +
      ( yy[,i,1] * ( 1 - zz[,2] ) ) +
      ( yy[,j,2] * ( 1 - zz[,1] ) ) * ( (eta[i]) / (colSums( nipij )[j]) ) +
      ( 1 - zz[,1] ) * ( 1 - zz[,2] ) * ( (eta[i]) / sum( nipij ) )
  }

  # lambda2 restriction
  lambda2 <-
    -( rowSums( Nij ) + Ri +
         rowSums( sweep( nipij , 2 , Cj / colSums( nipij ) , "*" ) ) +
         M * eta ) / N
  a.pij <- sweep( a.pij , 2 , lambda2 , "+" )

  # pij = zero
  pij.zero.mat <- which( pij == 0 , arr.ind = TRUE )
  for ( k in seq_len( nrow( pij.zero.mat ) ) ) a.pij[ , pij.zero.mat[k,1] , pij.zero.mat[k,2] ] <- 0

  # coerce to matrix
  u.pij <- matrix( 0 , nrow = dim( a.pij )[1] , ncol = K^2 , byrow = TRUE )
  for ( i in seq_len( nrow( Nij ) ) ) u.pij[,Kmat[i,]] <- a.pij[,i,]

  ### matrix of linearized variables

  # build Umat
  Umat <- do.call( cbind , list( u.psi , u.rho , u.tau , u.eta , u.pij ) )

  # test equality (within some tolerance)
  stopifnot( all.equal( colSums( Umat * ww ) , modelA.WVec( this.theta , Amat ) , check.attributes = FALSE , scale = Inf ) )

  ### calculate jacobian

  # jacobian matrix
  Jmat <- numDeriv::jacobian( modelA.WVec , this.theta , method = "complex" , side = NULL , CountMatrix = Amat )

  # inverse of the jacobian matrix
  Jmat.inv <- MASS::ginv( Jmat )

  # calculate variance
  Umat.adj <- t( apply( Umat , 1 , function(z) crossprod( t(Jmat.inv) , z ) ) )
  u.psi <- Umat.adj[ , 1 ]
  u.rho <- Umat.adj[ , 2 ]
  u.tau <- Umat.adj[ , 3 ]
  u.eta <- Umat.adj[ , 3 + seq_len(K) ]
  u.pij <- Umat.adj[ , (K+3) + seq_len(K^2) ]
  u.pij[ , which( t( pij ) == 0 , arr.ind = FALSE ) ] <- 0
  a.pij <- array( 0 , dim = c( nrow( xx ) , nrow( Nij ) , ncol( Nij ) ) )
  for ( i in seq_len( ncol( Nij ) ) ) {
    a.pij[,i,] <- u.pij[ , Kmat[ i , ] ]
  }

  ##### other variances

  # net flows
  u.nipij <- array( 0 , dim = c( nrow( xx ) , nrow( Nij ) , ncol( Nij ) ) )
  for ( i in seq_len( nrow(Nij) ) ) for ( j in seq_len( ncol( Nij ) ) ) {
    u.nipij[,i,j] <- ( pij[i,j] * u.eta[,i] + eta[i] * a.pij[,i,j] )
  }

  # gross flows
  a.muij <- array( 0 , dim = c( nrow( xx ) , nrow( Nij ) , ncol( Nij ) ) )
  for ( i in seq_len( nrow(Nij) ) ) for ( j in seq_len( ncol( Nij ) ) ) {
    a.muij[,i,j] <- N * u.nipij[,i,j] + nipij[i,j]
  }
  u.muij <- matrix( 0 , nrow = dim( a.muij )[1] , ncol = K^2 , byrow = TRUE )
  for ( i in seq_len( nrow( Nij ) ) ) u.muij[,Kmat[i,]] <- a.muij[,i,]
  rm( a.muij )

  # final distribution
  u.gamma <- apply( u.nipij , c(1,3) , sum )
  for ( j in seq_len( nrow( Nij ) ) ) u.gamma[,j] <- rowSums( u.nipij[,,j] )

  # delta
  delta <- N * ( colSums( nipij ) - eta )
  u.delta <- sweep( N * ( u.gamma - u.eta ) , 2 , ( colSums( nipij ) - eta ) , "+" )

  ##### split full matrix

  # build list of linearized variables
  llin <-
    list(
      "psi" = u.psi ,
      "rho" = u.rho ,
      "tau" = u.tau ,
      "eta" = u.eta ,
      "pij" = u.pij ,
      "muij" = u.muij ,
      "gamma" = u.gamma ,
      "delta" = u.delta )

  # return list
  return( llin )

}


