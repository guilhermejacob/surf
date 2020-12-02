### pseudo-likelihood functions

# Model A
pll.modA <- function( psi , rho , tau , eta , pij , CountMatrix ) {

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )

  # matrix blocks
  Part.Nij <- nipij * psi * rho
  Part.Cj <- colSums( nipij * ( 1 - psi ) * ( 1 - tau ) )
  Part.Ri <- rowSums( nipij * psi * ( 1 - rho ) )
  Part.M <- sum( nipij * ( 1 - psi ) * tau )

  # build matrix
  expected.props <- rbind( cbind( Part.Nij , Part.Cj ) , c( Part.Ri , Part.M ) )

  # evaluate
  sum( CountMatrix * log( expected.props ) ) # unconstrained

}

# Model B
pll.modB <- function( psi , rho , tau , eta , pij , CountMatrix ) {

  # intermediate computations
  nipij <- sweep( pij , 1 , eta , "*" )
  psinipij <- sweep( nipij , 1 , psi , "*" )
  psicnipij <- sweep( pij , 1 , 1 - psi , "*" )

  # matrix blocks
  Part.Nij <- psinipij * rho
  Part.Cj <- colSums( psicnipij * ( 1 - tau ) )
  Part.Ri <- rowSums( psinipij * ( 1 - rho ) )
  Part.M <- sum( psicnipij * tau )

  # build matrix
  expected.props <- rbind( cbind( Part.Nij , Part.Cj ) , c( Part.Ri , Part.M ) )

  # evaluate
  sum( CountMatrix * log( expected.props ) ) # unconstrained

}

# Model C
pll.modC <- function( psi , rho , tau , eta , pij , CountMatrix ) {

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

# Modelo D
pll.modD <- function( psi , rho , tau , eta , pij , CountMatrix ) {

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
  expected.props <- rbind( cbind( Part.Nij , Part.Cj ) , c( Part.Ri , Part.M ) )

  # evaluate
  sum( CountMatrix * log( expected.props ) ) # unconstrained

}
