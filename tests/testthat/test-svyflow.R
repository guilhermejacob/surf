context("svyflow estimates with nonresponse")

# load libraries
library( survey )
library( surf )

# define sample and population sizes
n <- as.integer( 10^4 )
N <- as.integer( 10^5 )

# initial and transition probabilities for parametric model
eta.pop <- c( .9 , .05, .05 )
pij.pop <- matrix( c(.80, .15, .05, .30, .60, .10, .10, .10, .80 ) , ncol = 3 , nrow = 3 , byrow = T )

# expected gross flows under parametric model
muij.pop <- N * sweep( pij.pop , 1 , eta.pop , "*" )

# parameters for non-reponse mechanism
psi.pop <- .80
rho.pop <- c( .80 , .50 , .50 )
tau.pop <- c( .90 , .70 , .60 )

# classification labels matrix
class_table <- expand.grid( data.frame( y0 = seq_len( nrow(pij.pop) ) , y1 = seq_len( nrow(pij.pop) ) ) )
class_table <- class_table[ order( class_table$y0 ) , ]
class_table[ ,"k_ij" ] <- as.character( seq_len( nrow( class_table ) ) )

### GENERATE POPULATION FROM PARAMETRIC MODEL

# generate N observations from parametric model
pop_fullresponse <- t( rmultinom( N , size = 1 , prob = as.numeric( t( sweep( pij.pop , 1 , eta.pop , "*" ) ) ) ) )

# transforms into transition typologies
pop_fullresponse <- apply( pop_fullresponse , 1 , function( z ) seq_len( ncol( pop_fullresponse ) )[ as.logical( z ) ] )

# transforms into categories
pop_fullresponse <- data.frame( "id" = seq_len( N ) , "k_ij" = pop_fullresponse , row.names = NULL , stringsAsFactors = FALSE )

# combine with classification
pop_fullresponse <- merge( pop_fullresponse , class_table , by = c( "k_ij" ) , all.x = TRUE , all.y = FALSE , sort = TRUE )

# reorder observations
pop_fullresponse <- pop_fullresponse[ order( pop_fullresponse$id ) , ]

# remove rownames
rownames( pop_fullresponse ) <- NULL

# transforms into factors
pop_fullresponse[, c( "y0" , "y1" ) ] <- lapply( pop_fullresponse[, c( "y0" , "y1" ) ] , factor , levels = c( 1:3 ) , labels = 1:3 )

### STOCHASTIC NON-RESPONSE MECHANISM

# copies full response population data
pop_nonrespose <- pop_fullresponse

# adds non-response in time t-1
pop_nonrespose[ as.logical( rbinom( nrow( pop_fullresponse ) , 1 , 1 - psi.pop ) ) , "y0" ] <- NA

# adds response transition in time t-1
for ( i in seq_along( levels( pop_fullresponse$y0 ) ) ) {
  pop_nonrespose[ pop_fullresponse$y0 == levels( pop_fullresponse$y0 )[i] &
                    !is.na( pop_nonrespose$y0 ) &
                    as.logical( rbinom( nrow( pop_fullresponse ) , 1 , 1 - rho.pop[i] ) ) , "y1" ] <- NA
}

# adds non-response transition in time t
for ( i in seq_along( levels( pop_fullresponse$y0 ) ) ) {
  pop_nonrespose[ pop_fullresponse$y0 == levels( pop_fullresponse$y0 )[i] &
                    is.na( pop_nonrespose$y0 ) &
                    as.logical( rbinom( nrow( pop_fullresponse ) , 1 , tau.pop[i] ) ) , "y1" ] <- NA
}

### SELECT SAMPLE

# select n observations from SRS
smp.df <- pop_nonrespose[ sample( N , n ) , c( "y0" , "y1" ) ]

# create sample design info
smp.df$prob <- n / N # selection probability
smp.df$fpcs <- N     # finite population correction

# build survey design object
flowdes.srs <-
  svydesign( ids = ~ 1 ,
             probs = ~ prob ,
             data = smp.df ,
             nest = TRUE )

# # create replicate designs
# flowdes_srs_rep <- as.surfrdesign( flowdes_srs , type = "bootstrap" , replicate = 50 )

# remove objects
rm( pop_fullresponse , pop_nonrespose , smp.df , class_table ) ; gc()

# estimate flows
flow_srs_lin <- svyflow( ~y0+y1 , flowdes.srs , model = "C" )
# flow_srs_rep <- svyflow( ~y , flowdes_srs_rep , model = "C" )

# test extraction of associated measures
test_that( "extraction of estimates" , {

  # point estimates
  expect_identical( coef( flow_srs_lin$psi ) , surf:::coef.svymstat( flow_srs_lin$psi ) )
  expect_identical( coef( flow_srs_lin$rho ) , surf:::coef.svymstat( flow_srs_lin$rho ) )
  expect_identical( coef( flow_srs_lin$tau ) , surf:::coef.svymstat( flow_srs_lin$tau ) )
  expect_identical( coef( flow_srs_lin$eta ) , surf:::coef.svymstat( flow_srs_lin$eta ) )
  expect_identical( coef( flow_srs_lin$muij ) , surf:::coef.svymstat( flow_srs_lin$muij ) )
  expect_identical( coef( flow_srs_lin$pij ) , surf:::coef.svymstat( flow_srs_lin$pij ) )
  # expect_identical( coef( flow_srs_rep$psi ) , surf:::coef.svymstat( flow_srs_rep$psi ) )
  # expect_identical( coef( flow_srs_rep$rho ) , surf:::coef.svymstat( flow_srs_rep$rho ) )
  # expect_identical( coef( flow_srs_rep$tau ) , surf:::coef.svymstat( flow_srs_rep$tau ) )
  # expect_identical( coef( flow_srs_rep$eta ) , surf:::coef.svymstat( flow_srs_rep$eta ) )
  # expect_identical( coef( flow_srs_rep$muij ) , surf:::coef.svymstat( flow_srs_rep$muij ) )
  # expect_identical( coef( flow_srs_rep$pij ) , surf:::coef.svymstat( flow_srs_rep$pij ) )

  # variances
  expect_identical( vcov( flow_srs_lin$psi ) , survey:::vcov.svystat( flow_srs_lin$psi ) )
  expect_identical( vcov( flow_srs_lin$rho ) , survey:::vcov.svystat( flow_srs_lin$rho ) )
  expect_identical( vcov( flow_srs_lin$tau ) , survey:::vcov.svystat( flow_srs_lin$tau ) )
  expect_identical( vcov( flow_srs_lin$eta ) , survey:::vcov.svystat( flow_srs_lin$eta ) )
  expect_identical( vcov( flow_srs_lin$muij ) , surf:::vcov.svymstat( flow_srs_lin$muij ) )
  expect_identical( vcov( flow_srs_lin$pij ) , surf:::vcov.svymstat( flow_srs_lin$pij ) )
  # expect_identical( vcov( flow_srs_rep$psi ) , survey:::vcov.svystat( flow_srs_rep$psi ) )
  # expect_identical( vcov( flow_srs_rep$rho ) , survey:::vcov.svystat( flow_srs_rep$rho ) )
  # expect_identical( vcov( flow_srs_rep$tau ) , survey:::vcov.svystat( flow_srs_rep$tau ) )
  # expect_identical( vcov( flow_srs_rep$eta ) , survey:::vcov.svystat( flow_srs_rep$eta ) )
  # expect_identical( vcov( flow_srs_rep$muij ) , surf:::vcov.svymstat( flow_srs_rep$muij ) )
  # expect_identical( vcov( flow_srs_rep$pij ) , surf:::vcov.svymstat( flow_srs_rep$pij ) )

  # standard errors
  expect_identical( SE( flow_srs_lin$psi ) , survey:::SE.svystat( flow_srs_lin$psi ) )
  expect_identical( SE( flow_srs_lin$rho ) , survey:::SE.svystat( flow_srs_lin$rho ) )
  expect_identical( SE( flow_srs_lin$tau ) , survey:::SE.svystat( flow_srs_lin$tau ) )
  expect_identical( SE( flow_srs_lin$eta ) , survey:::SE.svystat( flow_srs_lin$eta ) )
  expect_identical( SE( flow_srs_lin$muij ) , surf:::SE.svymstat( flow_srs_lin$muij ) )
  expect_identical( SE( flow_srs_lin$pij ) , surf:::SE.svymstat( flow_srs_lin$pij ) )
  # expect_identical( SE( flow_srs_rep$psi ) , survey:::SE.svystat( flow_srs_rep$psi ) )
  # expect_identical( SE( flow_srs_rep$rho ) , survey:::SE.svystat( flow_srs_rep$rho ) )
  # expect_identical( SE( flow_srs_rep$tau ) , survey:::SE.svystat( flow_srs_rep$tau ) )
  # expect_identical( SE( flow_srs_rep$eta ) , survey:::SE.svystat( flow_srs_rep$eta ) )
  # expect_identical( SE( flow_srs_rep$muij ) , surf:::SE.svymstat( flow_srs_rep$muij ) )
  # expect_identical( SE( flow_srs_rep$pij ) , surf:::SE.svymstat( flow_srs_rep$pij ) )

} )

# test against bias
test_that("compare point estimates vs population values",{

  # linearized design
  expect_equivalent( coef( flow_srs_lin$psi ) , psi.pop , tolerance = .50 )
  expect_equivalent( coef( flow_srs_lin$rho ) , rho.pop , tolerance = .50 )
  expect_equivalent( coef( flow_srs_lin$tau ) , tau.pop , tolerance = .50 )
  expect_equivalent( coef( flow_srs_lin$eta ) , eta.pop , tolerance = .20 )
  expect_equivalent( coef( flow_srs_lin$muij ) , muij.pop , tolerance = .20 )

  # # replicate design
  # expect_equivalent( coef( flow_srs_rep$psi ) , psi.pop , tolerance = .20 )
  # expect_equivalent( coef( flow_srs_rep$rho ) , rho.pop , tolerance = .20 )
  # expect_equivalent( coef( flow_srs_rep$tau ) , tau.pop , tolerance = .20 )
  # expect_equivalent( coef( flow_srs_rep$eta ) , eta.pop , tolerance = .20 )
  # expect_equivalent( coef( flow_srs_rep$muij ) , muij.pop , tolerance = .20 )

} )

# # linearized vs replicate
# test_that("compare linearized vs replicate: point estimates",{
#   expect_identical( coef( flow_srs_lin$psi ) , coef( flow_srs_rep$psi ) )
#   expect_identical( coef( flow_srs_lin$rho ) , coef( flow_srs_rep$rho ) )
#   expect_identical( coef( flow_srs_lin$tau ) , coef( flow_srs_rep$tau ) )
#   expect_identical( coef( flow_srs_lin$eta ) , coef( flow_srs_rep$eta ) )
#   expect_identical( coef( flow_srs_lin$muij ) , coef( flow_srs_rep$muij ) )
# } )
#
# test_that("compare linearized vs replicate: standard errors",{
#   expect_equivalent( SE( flow_srs_lin$psi ) , SE( flow_srs_rep$psi ) , tolerance = .05 )
#   expect_equivalent( SE( flow_srs_lin$rho ) , SE( flow_srs_rep$rho ) , tolerance = .05 )
#   expect_equivalent( SE( flow_srs_lin$tau ) , SE( flow_srs_rep$tau ) , tolerance = .05 )
#   expect_equivalent( SE( flow_srs_lin$eta ) , SE( flow_srs_rep$eta ) , tolerance = .05 )
#   expect_equivalent( SE( flow_srs_lin$muij ) , SE( flow_srs_rep$muij ) , tolerance = .30 )
# } )

