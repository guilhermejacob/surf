# function for model variance calculation
rao.scott <- function( xx , visit.number , ww , model.estimates , design , n.parms ) {

  # calculate degrees of freedom
  K <- dim( model.estimates )[1] - 1
  V <- dim( model.estimates )[3]
  n.counts <- sum( model.estimates > 0 )
  n.restr <- if ( V > 1 ) K+2 else K+1
  n.parms <- if ( V > 1 ) n.parms + V else n.parms
  df <- n.counts + n.restr - n.parms - 1

  # add non-reponse
  xx[,] <- lapply( xx , addNA )
  xx <- xx[, rev( seq_len( ncol( xx ) ) ) ]

  # create vector
  Avec <- interaction( do.call( cbind , list( xx , visit.number ) ) , drop = FALSE )
  Avec <- data.frame( FC = Avec )
  Amat <- stats::model.matrix( ~0+. , data = Avec[ ,, drop = FALSE ]  )

  # collect sample size
  smalln <- sum( ww > 0 )

  # calculate observed proportions
  observed.props <- svymean( Amat , design , na.rm = FALSE )

  # calculate chi distances
  est.vec <- c( aperm( model.estimates , c(2,1,3) ) )
  obs.vec <- as.numeric( coef( observed.props ) )
  chimat <- ( est.vec - obs.vec )^2 / est.vec
  obs.mat <- array( obs.vec , dim = dim( model.estimates ) )

  # # calculate deff
  # Vdes <- vcov( observed.props )
  # Dmat <- diag( obs.vec )
  # iDmat <- diag(ifelse(obs.vec == 0, 0, 1/obs.vec))
  # Vsrs <- ( Dmat - outer( obs.vec , obs.vec ) ) / (smalln-1)
  # # Delta <- Vdes / Vsrs
  # delta.value <- ( diag(Vdes) / diag(Vsrs) )[diag(Vsrs) > 0]

  # unadjusted chi-square
  statistic <- smalln * sum( chimat ) / df

  # calculate deff
  Vdes <- vcov( observed.props )
  Dmat <- diag( obs.vec )
  iDmat <- diag(ifelse(obs.vec == 0, 0, 1/obs.vec))
  Vsrs <- ( Dmat - outer( obs.vec , obs.vec ) ) / (smalln-1)
  Delta <- crossprod( MASS::ginv( Vsrs ) , Vdes )
  # delta.value <- as.numeric( eigen( Delta , only.values = TRUE )$values )
  delta.value <- diag( Delta )

  # first-order correction
  delta.bar <- sum( delta.value ) / df
  statistic <- statistic / delta.bar

  # second-order correction
  d0 <- sum(diag(Delta))^2/(sum(diag(Delta %*% Delta)))
  nu <- degf( design)

  # score
  warn <- options(warn = -1)
  pearson <- chisq.test( apply( obs.mat , 1:2 , sum ) )
  pearson$statistic <- statistic
  pearson$p.value <- pf( statistic , d0 , d0 * nu , lower.tail = FALSE )
  attr(pearson$statistic, "names") <- "F"
  pearson$parameter <- c(ndf = d0, ddf = d0 * nu)
  pearson$data.name <- "Expected vs. Observed Proportions"
  pearson$method <- "Pearson's X^2: Rao & Scott adjustment"

  # return
  pearson

}
