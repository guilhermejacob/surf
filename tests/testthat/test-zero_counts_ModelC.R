context("svyflow estimates with nonresponse and zero counts: Model C")

# set seed of random number generator
set.seed( 123 )

# load libraries
library( survey )
library( surf )
library( testthat )

# define population and sample size
N <- as.integer( 10^5 )
n <- as.integer( 10^4 )

# superpopulation hyperparameters
eta.pop <- c( .30 , .30, .20, .20 )
pij.pop <- matrix( c(.60,  .10 , .20, .10,
                     .30, .50, 0, .20,
                     .20, .20, .30, .30,
                     0 , .20 , .30, .50 ) , nrow = 4 , byrow = TRUE )

# expected gross flows
muij.pop <- N * sweep( pij.pop , 1 , eta.pop , "*" )

# non-response pattern
psi.pop <- .70
rho.pop <- c( .80 , .50 , .70 , .60 )
tau.pop <- c( .80 , .70 , .60 , .80 )

# population parameters
poplist <-
  list( "psi" = psi.pop ,
        "rho" = rho.pop ,
        "tau" =  tau.pop ,
        "eta" = eta.pop ,
        "gamma" = colSums( sweep( pij.pop , 1 , eta.pop , "*" ) ) ,
        "pij" = pij.pop ,
        "muij" = muij.pop )

# classification matrix
state.table <- expand.grid( data.frame( v0 = seq_len( nrow(pij.pop) +1 ) , v1 = seq_len( nrow(pij.pop) + 1 ) ) )
state.table <- state.table[ order( state.table$v0 ) , ]
state.table[ ,"k_ij" ] <- as.character( seq_len( nrow( state.table ) ) )
state.table$v0[ state.table$v0 == 5 ] <- NA
state.table$v1[ state.table$v1 == 5 ] <- NA

# intermediate computations
nipij <- sweep( pij.pop , 1 , eta.pop , "*" )
rhonipij <- sweep( nipij , 1 , rho.pop , "*" )
rhocnipij <- sweep( nipij , 1 , 1 - rho.pop , "*" )
taunipij <- sweep( nipij , 1 , tau.pop , "*" )
taucnipij <- sweep( nipij , 1 , 1 - tau.pop , "*" )

# matrix blocks
Part.Nij <- psi.pop * rhonipij
Part.Ri <- rowSums( psi.pop * rhocnipij )
Part.Cj <- colSums( (1 - psi.pop) * taucnipij )
Part.M <- sum( (1 - psi.pop ) * taunipij )

# build matrix
exp.props <- rbind( cbind( Part.Nij , Part.Ri ) , c( Part.Cj , Part.M ) )
dimnames( exp.props ) <- NULL

# extract sample
smp.df <- t( rmultinom( n , size = 1 , prob = as.numeric( t( exp.props ) ) ) )
smp.df <- apply( smp.df , 1 , function( z ) seq_len( ncol( smp.df ) )[ as.logical( z ) ] )
smp.df <- data.frame( "id" = seq_len( n ) , "k_ij" = smp.df , row.names = NULL , stringsAsFactors = FALSE )
smp.df <- merge( smp.df , state.table , by = c( "k_ij" ) , all.x = TRUE , all.y = FALSE , sort = FALSE )
smp.df <- smp.df[ order( smp.df$id ) , ]
rownames( smp.df ) <- NULL
smp.df[, c( "v0" , "v1" ) ] <- lapply( smp.df[, c( "v0" , "v1" ) ] , factor , levels = c( 1:4 ) , labels = 1:4 )
table( smp.df[ , c("v0","v1")] , useNA = "always" )

# sampling des.lin information
smp.df$wgt    <- N / n
smp.df$fpcs   <- N

# declare survey des.lin object
des.lin <-
  svydesign( ids = ~ 1 ,
             weights = ~ wgt ,
             fpc = ~fpcs ,
             data = smp.df ,
             nest = TRUE )

# resampling object
des.rep <- as.svrepdesign( des.lin , "bootstrap" , replicates = 100 )

# estima contagens
svytable( ~v0+v1 , des.lin , addNA = TRUE )

# estimate gross flows
suppressWarnings( flow_srs_lin <- svyflow( ~v0+v1 , des.lin , model = "C" , verbose = FALSE , as.zero.flows = TRUE ) )
suppressWarnings( flow_srs_rep <- svyflow( ~v0+v1 , des.rep , model = "C" , verbose = FALSE , as.zero.flows = TRUE ) )

# test extraction of associated measures
test_that( "extraction of estimates: linearization" , {

  # point estimates
  expect_identical( coef( flow_srs_lin$psi ) , survey:::coef.svystat( flow_srs_lin$psi ) )
  expect_identical( coef( flow_srs_lin$rho ) , survey:::coef.svystat( flow_srs_lin$rho ) )
  expect_identical( coef( flow_srs_lin$tau ) , survey:::coef.svystat( flow_srs_lin$tau ) )
  expect_identical( coef( flow_srs_lin$eta ) , survey:::coef.svystat( flow_srs_lin$eta ) )
  expect_identical( coef( flow_srs_lin$pij ) , surf:::coef.svymstat( flow_srs_lin$pij ) )
  expect_identical( coef( flow_srs_lin$muij ) , surf:::coef.svymstat( flow_srs_lin$muij ) )

  # variances
  expect_identical( vcov( flow_srs_lin$psi ) , survey:::vcov.svystat( flow_srs_lin$psi ) )
  expect_identical( vcov( flow_srs_lin$rho ) , survey:::vcov.svystat( flow_srs_lin$rho ) )
  expect_identical( vcov( flow_srs_lin$tau ) , survey:::vcov.svystat( flow_srs_lin$tau ) )
  expect_identical( vcov( flow_srs_lin$eta ) , survey:::vcov.svystat( flow_srs_lin$eta ) )
  expect_identical( vcov( flow_srs_lin$pij ) , surf:::vcov.svymstat( flow_srs_lin$pij ) )
  expect_identical( vcov( flow_srs_lin$muij ) , surf:::vcov.svymstat( flow_srs_lin$muij ) )

  # standard errors
  expect_identical( SE( flow_srs_lin$psi ) , survey:::SE.svystat( flow_srs_lin$psi ) )
  expect_identical( SE( flow_srs_lin$rho ) , survey:::SE.svystat( flow_srs_lin$rho ) )
  expect_identical( SE( flow_srs_lin$tau ) , survey:::SE.svystat( flow_srs_lin$tau ) )
  expect_identical( SE( flow_srs_lin$eta ) , survey:::SE.svystat( flow_srs_lin$eta ) )
  expect_identical( SE( flow_srs_lin$pij ) , surf:::SE.svymstat( flow_srs_lin$pij ) )
  expect_identical( SE( flow_srs_lin$muij ) , surf:::SE.svymstat( flow_srs_lin$muij ) )

} )

test_that( "extraction of estimates: resampling" , {

  # point estimates
  expect_identical( coef( flow_srs_rep$psi ) , survey:::coef.svystat( flow_srs_rep$psi ) )
  expect_identical( coef( flow_srs_rep$rho ) , survey:::coef.svystat( flow_srs_rep$rho ) )
  expect_identical( coef( flow_srs_rep$tau ) , survey:::coef.svystat( flow_srs_rep$tau ) )
  expect_identical( coef( flow_srs_rep$eta ) , survey:::coef.svystat( flow_srs_rep$eta ) )
  expect_identical( coef( flow_srs_rep$pij ) , surf:::coef.svymstat( flow_srs_rep$pij ) )
  expect_identical( coef( flow_srs_rep$muij ) , surf:::coef.svymstat( flow_srs_rep$muij ) )

  # variances
  expect_identical( vcov( flow_srs_rep$psi ) , survey:::vcov.svystat( flow_srs_rep$psi ) )
  expect_identical( vcov( flow_srs_rep$rho ) , survey:::vcov.svystat( flow_srs_rep$rho ) )
  expect_identical( vcov( flow_srs_rep$tau ) , survey:::vcov.svystat( flow_srs_rep$tau ) )
  expect_identical( vcov( flow_srs_rep$eta ) , survey:::vcov.svystat( flow_srs_rep$eta ) )
  expect_identical( vcov( flow_srs_rep$pij ) , surf:::vcov.svymstat( flow_srs_rep$pij ) )
  expect_identical( vcov( flow_srs_rep$muij ) , surf:::vcov.svymstat( flow_srs_rep$muij ) )

  # standard errors
  expect_identical( SE( flow_srs_rep$psi ) , survey:::SE.svystat( flow_srs_rep$psi ) )
  expect_identical( SE( flow_srs_rep$rho ) , survey:::SE.svystat( flow_srs_rep$rho ) )
  expect_identical( SE( flow_srs_rep$tau ) , survey:::SE.svystat( flow_srs_rep$tau ) )
  expect_identical( SE( flow_srs_rep$eta ) , survey:::SE.svystat( flow_srs_rep$eta ) )
  expect_identical( SE( flow_srs_rep$pij ) , surf:::SE.svymstat( flow_srs_rep$pij ) )
  expect_identical( SE( flow_srs_rep$muij ) , surf:::SE.svymstat( flow_srs_rep$muij ) )

} )

# test against bias
test_that("compare point estimates vs population values",{

  # linearized des.lin
  expect_equivalent( coef( flow_srs_lin$psi ) , psi.pop , tolerance = .30 )
  expect_equivalent( coef( flow_srs_lin$rho ) , rho.pop , tolerance = .30 )
  expect_equivalent( coef( flow_srs_lin$tau ) , tau.pop , tolerance = .30 )
  expect_equivalent( coef( flow_srs_lin$eta ) , eta.pop , tolerance = .10 )
  expect_equivalent( coef( flow_srs_lin$pij , to.matrix = TRUE ) , pij.pop , tolerance = .10 )
  expect_equivalent( coef( flow_srs_lin$muij , to.matrix = TRUE ) , muij.pop , tolerance = .20 )

  # resampling des.lin
  expect_equivalent( coef( flow_srs_rep$psi ) , psi.pop , tolerance = .30 )
  expect_equivalent( coef( flow_srs_rep$rho ) , rho.pop , tolerance = .30 )
  expect_equivalent( coef( flow_srs_rep$tau ) , tau.pop , tolerance = .30 )
  expect_equivalent( coef( flow_srs_rep$eta ) , eta.pop , tolerance = .10 )
  expect_equivalent( coef( flow_srs_rep$pij , to.matrix = TRUE ) , pij.pop , tolerance = .10 )
  expect_equivalent( coef( flow_srs_rep$muij , to.matrix = TRUE ) , muij.pop , tolerance = .20 )

} )
