##### Model A #####

# pseudo log-likelihood
modelA.loglik <- function( theta , CountMatrix ) {

  # nummber of categories
  K <- dim( CountMatrix )[1] - 1

  # rebuild parameters form vector theta
  psi <- theta[ 1 ]
  rho <- theta[ 2 ]
  tau <- theta[ 3 ]
  eta <- theta[ 3 + seq_len( K ) ]
  pij <- theta[ seq( 4 + K , ( 4 + K ) + K^2 - 1 ) ]

  # rebuild matrices
  pij <- matrix( pij , ncol = K , nrow = K , byrow = TRUE )

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )

  # matrix blocks
  Part.Nij <- nipij * psi * rho
  Part.Cj <- colSums( nipij * ( 1 - psi ) * ( 1 - tau ) )
  Part.Ri <- rowSums( nipij * psi * ( 1 - rho ) )
  Part.M <- sum( nipij * ( 1 - psi ) * tau )

  # build matrix
  expected.props <- rbind( cbind( Part.Nij , Part.Ri ) , c( Part.Cj , Part.M ) )

  # evaluate
  sum( CountMatrix * log( expected.props ) ) # unconstrained

}


# expected proportions
modelA.expected <- function( theta , K ) {

  # rebuild parameters form vector theta
  psi <- theta[ 1 ]
  rho <- theta[ 2 ]
  tau <- theta[ 3 ]
  eta <- theta[ 3 + seq_len( K ) ]
  pij <- theta[ seq( 4 + K , ( 4 + K ) + K^2 - 1 ) ]

  # rebuild matrices
  pij <- matrix( pij , ncol = K , nrow = K , byrow = TRUE )

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )

  # matrix blocks
  Part.Nij <- nipij * psi * rho
  Part.Cj <- colSums( nipij * ( 1 - psi ) * ( 1 - tau ) )
  Part.Ri <- rowSums( nipij * psi * ( 1 - rho ) )
  Part.M <- sum( nipij * ( 1 - psi ) * tau )

  # build matrix
  return( rbind( cbind( Part.Nij , Part.Ri ) , c( Part.Cj , Part.M ) ) )

}

##### Model B #####

# pseudo log-likelihood
modelB.loglik <- function( theta , CountMatrix ) {

  # nummber of categories
  K <- dim( CountMatrix )[1] - 1

  # Obtain sample estimates of population flows as described in table 3.1 of Rojas et al. (2014)
  Nij <- as.matrix( CountMatrix[ -nrow( CountMatrix ) , -ncol( CountMatrix ) ] )
  Ri <- CountMatrix[ , ncol(CountMatrix) ][ - nrow( CountMatrix ) ]
  Cj <- CountMatrix[ nrow(CountMatrix) , ][ - ncol( CountMatrix ) ]
  M <- CountMatrix[ nrow( CountMatrix ) , ncol( CountMatrix ) ]
  N  <- sum( Nij ) + sum( Ri ) + sum( Cj ) + M

  # rebuild parameters form vector theta
  psi <- theta[ 1:K ]
  rho <- theta[ K+1 ]
  tau <- theta[ K+2 ]
  eta <- theta[ seq( K+3 , K+3 + ( K - 1 ) ) ]
  pij <- theta[ seq( 2*K+3 , ( 2*K+3 ) + ( K^2 - 1 ) ) ]

  # rebuild matrices
  pij <- matrix( pij , ncol = K , nrow = K , byrow = TRUE )

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )
  psinipij <- sweep( nipij , 1 , psi , "*" )
  psicnipij <- sweep( nipij , 1 , 1 - psi , "*" )

  # matrix blocks
  Part.Nij <- psinipij * rho
  Part.Cj <- colSums( psicnipij * ( 1 - tau ) )
  Part.Ri <- rowSums( psinipij * ( 1 - rho ) )
  Part.M <- sum( psicnipij * tau )

  # build matrix
  expected.props <- rbind( cbind( Part.Nij , Part.Ri ) , c( Part.Cj , Part.M ) )
  expected.props <- expected.props / sum( expected.props )

  # evaluate
  sum( CountMatrix * log( expected.props ) ) # unconstrained

}

# expected proportions
modelB.expected <- function( theta , K ) {

  # rebuild parameters form vector theta
  psi <- theta[ 1:K ]
  rho <- theta[ K+1 ]
  tau <- theta[ K+2 ]
  eta <- theta[ seq( K+3 , K+3 + ( K - 1 ) ) ]
  pij <- theta[ seq( 2*K+3 , ( 2*K+3 ) + ( K^2 - 1 ) ) ]

  # rebuild matrices
  pij <- matrix( pij , ncol = K , nrow = K , byrow = TRUE )

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )
  psinipij <- sweep( nipij , 1 , psi , "*" )
  psicnipij <- sweep( nipij , 1 , 1 - psi , "*" )

  # matrix blocks
  Part.Nij <- psinipij * rho
  Part.Cj <- colSums( psicnipij * ( 1 - tau ) )
  Part.Ri <- rowSums( psinipij * ( 1 - rho ) )
  Part.M <- sum( psicnipij * tau )

  # build matrix
  rbind( cbind( Part.Nij , Part.Ri ) , c( Part.Cj , Part.M ) )

}

##### Model C #####

# pseudo log-likelihood
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
  expected.props <- rbind( cbind( Part.Nij , Part.Ri ) , c( Part.Cj , Part.M ) )

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
  rbind( cbind( Part.Nij , Part.Ri ) , c( Part.Cj , Part.M ) )

}

##### Model D #####

# pseudo log-likelihood
modelD.loglik <- function( theta , CountMatrix ) {

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
  rhonipij <- sweep( nipij , 2 , rho , "*" )
  rhocnipij <- sweep( nipij , 2 , 1 - rho , "*" )
  taunipij <- sweep( nipij , 2 , tau , "*" )
  taucnipij <- sweep( nipij , 2 , 1 - tau , "*" )

  # matrix blocks
  Part.Nij <- psi * rhonipij
  Part.Cj <- colSums( (1 - psi) * taucnipij )
  Part.Ri <- rowSums( psi * rhocnipij )
  Part.M <- sum( (1 - psi ) * taunipij )

  # build matrix
  expected.props <- rbind( cbind( Part.Nij , Part.Ri ) , c( Part.Cj , Part.M ) )

  # evaluate
  sum( CountMatrix * log( expected.props ) ) # unconstrained

}


# expected.proportions
modelD.expected <- function( theta , K ) {

  # rebuild parameters form vector theta
  psi <- theta[ 1 ]
  rho <- theta[ 1+seq_len(K) ]
  tau <- theta[ (K+1) + seq_len(K) ]
  eta <- theta[ (2*K+1) + seq_len(K) ]
  pij <- theta[ 3*K+1 + seq_len(K^2) ]

  # test formats
  stopifnot( (3*K+1 + K^2 ) == length( theta ) )

  # rebuild matrices
  pij <- matrix( pij , ncol = K , nrow = K , byrow = TRUE )

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )
  rhonipij <- sweep( nipij , 2 , rho , "*" )
  rhocnipij <- sweep( nipij , 2 , 1 - rho , "*" )
  taunipij <- sweep( nipij , 2 , tau , "*" )
  taucnipij <- sweep( nipij , 2 , 1 - tau , "*" )

  # matrix blocks
  Part.Nij <- psi * rhonipij
  Part.Cj <- colSums( (1 - psi) * taucnipij )
  Part.Ri <- rowSums( psi * rhocnipij )
  Part.M <- sum( (1 - psi ) * taunipij )

  # build matrix
  rbind( cbind( Part.Nij , Part.Ri ) , c( Part.Cj , Part.M ) )

}
