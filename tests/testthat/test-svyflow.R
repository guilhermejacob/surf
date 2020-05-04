context("svyflow estimates with nonresponse")

# load libraries
library( surf )

# define sample and population sizes
n <- as.integer( 10^4 )
N <- as.integer( 10^5 )

# initial and transition probabilities for parametric model
eta_pop <- c( .9 , .05, .05 )
pij_pop <- matrix( c(.80, .15, .05, .30, .60, .10, .10, .10, .80 ) , ncol = 3 , nrow = 3 , byrow = T )

# expected gross flows under parametric model
muij_pop <- N * sweep( pij_pop , 1 , eta_pop , "*" )

# parameters for non-reponse mechanism
psi_pop <- .80
rhoRR_pop <- c( .80 , .50 , .50 )
rhoMM_pop <- c( .90 , .70 , .60 )

# classification labels matrix
class_table <- expand.grid( data.frame( v0 = seq_len( nrow(pij_pop) ) , v1 = seq_len( nrow(pij_pop) ) ) )
class_table <- class_table[ order( class_table$v0 ) , ]
class_table[ ,"k_ij" ] <- as.character( seq_len( nrow( class_table ) ) )

### GENERATE POPULATION FROM PARAMETRIC MODEL

# generate N observations from parametric model
pop_fullresponse <- t( rmultinom( N , size = 1 , prob = as.numeric( t( sweep( pij_pop , 1 , eta_pop , "*" ) ) ) ) )

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
pop_fullresponse[, c( "v0" , "v1" ) ] <- lapply( pop_fullresponse[, c( "v0" , "v1" ) ] , factor , levels = c( 1:3 ) , labels = 1:3 )

### STOCHASTIC NON-RESPONSE MECHANISM

# copies full response population data
pop_nonrespose <- pop_fullresponse

# adds non-response in time t-1
pop_nonrespose[ as.logical( rbinom( nrow( pop_fullresponse ) , 1 , 1 - psi_pop ) ) , "v0" ] <- NA

# adds response transition in time t-1
for ( i in seq_along( levels( pop_fullresponse$v0 ) ) ) {
  pop_nonrespose[ pop_fullresponse$v0 == levels( pop_fullresponse$v0 )[i] &
                    !is.na( pop_nonrespose$v0 ) &
                    as.logical( rbinom( nrow( pop_fullresponse ) , 1 , 1 - rhoRR_pop[i] ) ) , "v1" ] <- NA
}

# adds non-response transition in time t
for ( i in seq_along( levels( pop_fullresponse$v0 ) ) ) {
  pop_nonrespose[ pop_fullresponse$v0 == levels( pop_fullresponse$v0 )[i] &
                    is.na( pop_nonrespose$v0 ) &
                    as.logical( rbinom( nrow( pop_fullresponse ) , 1 , rhoMM_pop[i] ) ) , "v1" ] <- NA
}

### SELECT SAMPLE

# select n observations from SRS
smp_df <- pop_nonrespose[ sample( N , n ) , c( "v0" , "v1" ) ]

# create sample design info
smp_df$prob <- n / N            # selection probability
smp_df$fpcs <- N  # finite population correction

# build sample data frames
df0 <- smp_df[ , -2 , drop = FALSE ]
df1 <- smp_df[ , 2 , drop = FALSE ]
colnames( df0 )[1] <- "y"
colnames( df1 )[1] <- "y"

# build survey design object
flowdes_srs <-
  sfydesign( ids = ~ 1 ,
             probs = ~ prob ,
             data = list( df0 , df1 ) ,
             nest = TRUE )

# # create replicate designs
# flowdes_srs_rep <- as.surfrdesign( flowdes_srs , type = "bootstrap" , replicate = 50 )

# remove objects
rm( pop_fullresponse , pop_nonrespose , smp_df , df0 , df1 , class_table ) ; gc()

# estimate flows
flow_srs_lin <- svyflow( ~y , flowdes_srs , model = "C" )
# flow_srs_rep <- svyflow( ~y , flowdes_srs_rep , model = "C" )

# test extraction of associated measures
test_that( "extraction of estimates" , {

  # point estimates
  expect_identical( coef( flow_srs_lin$psi ) , surf:::coef.svymstat( flow_srs_lin$psi ) )
  expect_identical( coef( flow_srs_lin$rhoRR ) , surf:::coef.svymstat( flow_srs_lin$rhoRR ) )
  expect_identical( coef( flow_srs_lin$rhoMM ) , surf:::coef.svymstat( flow_srs_lin$rhoMM ) )
  expect_identical( coef( flow_srs_lin$eta ) , surf:::coef.svymstat( flow_srs_lin$eta ) )
  expect_identical( coef( flow_srs_lin$muij ) , surf:::coef.svymstat( flow_srs_lin$muij ) )
  expect_identical( coef( flow_srs_lin$pij ) , surf:::coef.svymstat( flow_srs_lin$pij ) )
  # expect_identical( coef( flow_srs_rep$psi ) , surf:::coef.svymstat( flow_srs_rep$psi ) )
  # expect_identical( coef( flow_srs_rep$rhoRR ) , surf:::coef.svymstat( flow_srs_rep$rhoRR ) )
  # expect_identical( coef( flow_srs_rep$rhoMM ) , surf:::coef.svymstat( flow_srs_rep$rhoMM ) )
  # expect_identical( coef( flow_srs_rep$eta ) , surf:::coef.svymstat( flow_srs_rep$eta ) )
  # expect_identical( coef( flow_srs_rep$muij ) , surf:::coef.svymstat( flow_srs_rep$muij ) )
  # expect_identical( coef( flow_srs_rep$pij ) , surf:::coef.svymstat( flow_srs_rep$pij ) )

  # variances
  expect_identical( vcov( flow_srs_lin$psi ) , survey:::vcov.svystat( flow_srs_lin$psi ) )
  expect_identical( vcov( flow_srs_lin$rhoRR ) , survey:::vcov.svystat( flow_srs_lin$rhoRR ) )
  expect_identical( vcov( flow_srs_lin$rhoMM ) , survey:::vcov.svystat( flow_srs_lin$rhoMM ) )
  expect_identical( vcov( flow_srs_lin$eta ) , survey:::vcov.svystat( flow_srs_lin$eta ) )
  expect_identical( vcov( flow_srs_lin$muij ) , surf:::vcov.svymstat( flow_srs_lin$muij ) )
  expect_identical( vcov( flow_srs_lin$pij ) , surf:::vcov.svymstat( flow_srs_lin$pij ) )
  # expect_identical( vcov( flow_srs_rep$psi ) , survey:::vcov.svystat( flow_srs_rep$psi ) )
  # expect_identical( vcov( flow_srs_rep$rhoRR ) , survey:::vcov.svystat( flow_srs_rep$rhoRR ) )
  # expect_identical( vcov( flow_srs_rep$rhoMM ) , survey:::vcov.svystat( flow_srs_rep$rhoMM ) )
  # expect_identical( vcov( flow_srs_rep$eta ) , survey:::vcov.svystat( flow_srs_rep$eta ) )
  # expect_identical( vcov( flow_srs_rep$muij ) , surf:::vcov.svymstat( flow_srs_rep$muij ) )
  # expect_identical( vcov( flow_srs_rep$pij ) , surf:::vcov.svymstat( flow_srs_rep$pij ) )

  # standard errors
  expect_identical( SE( flow_srs_lin$psi ) , survey:::SE.svystat( flow_srs_lin$psi ) )
  expect_identical( SE( flow_srs_lin$rhoRR ) , survey:::SE.svystat( flow_srs_lin$rhoRR ) )
  expect_identical( SE( flow_srs_lin$rhoMM ) , survey:::SE.svystat( flow_srs_lin$rhoMM ) )
  expect_identical( SE( flow_srs_lin$eta ) , survey:::SE.svystat( flow_srs_lin$eta ) )
  expect_identical( SE( flow_srs_lin$muij ) , surf:::SE.svymstat( flow_srs_lin$muij ) )
  expect_identical( SE( flow_srs_lin$pij ) , surf:::SE.svymstat( flow_srs_lin$pij ) )
  # expect_identical( SE( flow_srs_rep$psi ) , survey:::SE.svystat( flow_srs_rep$psi ) )
  # expect_identical( SE( flow_srs_rep$rhoRR ) , survey:::SE.svystat( flow_srs_rep$rhoRR ) )
  # expect_identical( SE( flow_srs_rep$rhoMM ) , survey:::SE.svystat( flow_srs_rep$rhoMM ) )
  # expect_identical( SE( flow_srs_rep$eta ) , survey:::SE.svystat( flow_srs_rep$eta ) )
  # expect_identical( SE( flow_srs_rep$muij ) , surf:::SE.svymstat( flow_srs_rep$muij ) )
  # expect_identical( SE( flow_srs_rep$pij ) , surf:::SE.svymstat( flow_srs_rep$pij ) )

} )

# test against bias
test_that("compare point estimates vs population values",{

  # linearized design
  expect_equivalent( coef( flow_srs_lin$psi ) , psi_pop , tolerance = .50 )
  expect_equivalent( coef( flow_srs_lin$rhoRR ) , rhoRR_pop , tolerance = .50 )
  expect_equivalent( coef( flow_srs_lin$rhoMM ) , rhoMM_pop , tolerance = .50 )
  expect_equivalent( coef( flow_srs_lin$eta ) , eta_pop , tolerance = .20 )
  expect_equivalent( coef( flow_srs_lin$muij ) , muij_pop , tolerance = .20 )

  # # replicate design
  # expect_equivalent( coef( flow_srs_rep$psi ) , psi_pop , tolerance = .20 )
  # expect_equivalent( coef( flow_srs_rep$rhoRR ) , rhoRR_pop , tolerance = .20 )
  # expect_equivalent( coef( flow_srs_rep$rhoMM ) , rhoMM_pop , tolerance = .20 )
  # expect_equivalent( coef( flow_srs_rep$eta ) , eta_pop , tolerance = .20 )
  # expect_equivalent( coef( flow_srs_rep$muij ) , muij_pop , tolerance = .20 )

} )

# # linearized vs replicate
# test_that("compare linearized vs replicate: point estimates",{
#   expect_identical( coef( flow_srs_lin$psi ) , coef( flow_srs_rep$psi ) )
#   expect_identical( coef( flow_srs_lin$rhoRR ) , coef( flow_srs_rep$rhoRR ) )
#   expect_identical( coef( flow_srs_lin$rhoMM ) , coef( flow_srs_rep$rhoMM ) )
#   expect_identical( coef( flow_srs_lin$eta ) , coef( flow_srs_rep$eta ) )
#   expect_identical( coef( flow_srs_lin$muij ) , coef( flow_srs_rep$muij ) )
# } )
#
# test_that("compare linearized vs replicate: standard errors",{
#   expect_equivalent( SE( flow_srs_lin$psi ) , SE( flow_srs_rep$psi ) , tolerance = .05 )
#   expect_equivalent( SE( flow_srs_lin$rhoRR ) , SE( flow_srs_rep$rhoRR ) , tolerance = .05 )
#   expect_equivalent( SE( flow_srs_lin$rhoMM ) , SE( flow_srs_rep$rhoMM ) , tolerance = .05 )
#   expect_equivalent( SE( flow_srs_lin$eta ) , SE( flow_srs_rep$eta ) , tolerance = .05 )
#   expect_equivalent( SE( flow_srs_lin$muij ) , SE( flow_srs_rep$muij ) , tolerance = .30 )
# } )

