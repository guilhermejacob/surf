context("subsets on survflow.design")

# random number generator seed
set.seed( 123 )

# population size
N = as.integer( 10^5 )

# sample size
n = as.integer(10^4)

# superpopulation model parameters
eta_pop <- c( .9 , .05, .05 )
pij_pop <- matrix( c(.80, .15, .05, .30, .60, .10, .10, .10, .80 ) , ncol = 3 , nrow = 3 , byrow = T )

# non-response parameters
psi_pop <- .8
rhoRR_pop <- .9
rhoMM_pop <- .7

# population gross flows
muij_pop <- N * sweep( pij_pop , 1 , eta_pop , "*" )

# classification matrix
class_table <- expand.grid( data.frame( v0 = seq_len( nrow(pij_pop) ) , v1 = seq_len( nrow(pij_pop) ) ) )
class_table <- class_table[ order( class_table$v0 ) , ]
class_table[ ,"k_ij" ] <- as.character( seq_len( nrow( class_table ) ) )

##### generate population from super-population model

# extract N-sized population
pop_fullresponse <- t( rmultinom( N , size = 1 , prob = as.numeric( t( sweep( pij_pop , 1 , eta_pop , "*" ) ) ) ) )

# apply transitions
pop_fullresponse <- apply( pop_fullresponse , 1 , function( z ) seq_len( ncol( pop_fullresponse ) )[ as.logical( z ) ] )
pop_fullresponse <- data.frame( "id" = seq_len( N ) , "k_ij" = pop_fullresponse , row.names = NULL , stringsAsFactors = FALSE )
pop_fullresponse <- merge( pop_fullresponse , class_table , by = c( "k_ij" ) , all.x = TRUE , all.y = FALSE , sort = TRUE )
pop_fullresponse <- pop_fullresponse[ order( pop_fullresponse$id ) , ]

# remove row names
rownames( pop_fullresponse ) <- NULL

# to factors
pop_fullresponse[, c( "v0" , "v1" ) ] <- lapply( pop_fullresponse[, c( "v0" , "v1" ) ] , factor , levels = c( 1:3 ) , labels = 1:3 )

##### apply non-response

# copy full resopnse population
pop_nonrespose <- pop_fullresponse

# adds non-response in time t-1
pop_nonrespose[ as.logical( rbinom( nrow( pop_fullresponse ) , 1 , 1 - psi_pop ) ) , "v0" ] <- NA

# adds non-response in time t for respondents in time t-1
pop_nonrespose[ !is.na( pop_nonrespose[,"v0"] ) & as.logical( rbinom( nrow( pop_nonrespose ) , 1 , 1 - rhoRR_pop ) ) , "v1" ] <- NA

# adds non-response in time t for non-respondents in time t-1
pop_nonrespose[  is.na( pop_nonrespose[,"v0"] ) & as.logical( rbinom( nrow( pop_nonrespose ) , 1, rhoMM_pop ) ) , "v1" ] <- NA

### extract sample from population

# extract SRS of size n
smp_df <- pop_nonrespose[ sample( N , n ) , c( "v0" , "v1" ) ]

# create sampling info
smp_df$prob <- n / N # selection probabilities

# adjust data.frames
df0 <- smp_df[ , -2 , drop = FALSE ]
df1 <- smp_df[ , 2 , drop = FALSE ]
colnames( df0 )[1] <- "y"
colnames( df1 )[1] <- "y"

# cria variável
df0$v1 <- rbinom( n , 1 , .5 ) + 1

# build survey design object
flowdes_srs <- sfydesign( ids = ~ 1 ,
                          probs = ~ prob ,
                          data = list( df0 , df1 ) ,
                          nest = TRUE )

# create replicate designs
flowdes_srs_rep <- as.surfrdesign( flowdes_srs , type = "bootstrap" , replicate = 50 )

# remove objects
rm( pop_fullresponse , pop_nonrespose , smp_df , df0 , df1 , class_table ) ; gc()

# estimate flows
flow_srs_lin <- svyflow( ~y , flowdes_srs )
flow_lin_sub1 <- svyflow( ~y , subset( flowdes_srs , v1 == 1 ) , na.rm = TRUE )
flow_lin_sub2 <- svyflow( ~y , subset( flowdes_srs , v1 == 2 ) , na.rm = TRUE )
flow_srs_rep <- svyflow( ~y , flowdes_srs_rep )
flow_rep_sub1 <- svyflow( ~y , subset( flowdes_srs_rep , v1 == 1 ) , na.rm = TRUE )
flow_rep_sub2 <- svyflow( ~y , subset( flowdes_srs_rep , v1 == 2 ) , na.rm = TRUE )

# test extraction of associated measures
test_that( "extraction of estimates" , {

  # point estimates
  expect_identical( coef( flow_lin_sub1$psi ) , surf:::coef.svymstat( flow_lin_sub1$psi ) )
  expect_identical( coef( flow_lin_sub1$rhoRR ) , surf:::coef.svymstat( flow_lin_sub1$rhoRR ) )
  expect_identical( coef( flow_lin_sub1$rhoMM ) , surf:::coef.svymstat( flow_lin_sub1$rhoMM ) )
  expect_identical( coef( flow_lin_sub1$eta ) , surf:::coef.svymstat( flow_lin_sub1$eta ) )
  expect_identical( coef( flow_lin_sub1$muij ) , surf:::coef.svymstat( flow_lin_sub1$muij ) )
  expect_identical( coef( flow_lin_sub1$pij ) , surf:::coef.svymstat( flow_lin_sub1$pij ) )
  expect_identical( coef( flow_lin_sub2$psi ) , surf:::coef.svymstat( flow_lin_sub2$psi ) )
  expect_identical( coef( flow_lin_sub2$rhoRR ) , surf:::coef.svymstat( flow_lin_sub2$rhoRR ) )
  expect_identical( coef( flow_lin_sub2$rhoMM ) , surf:::coef.svymstat( flow_lin_sub2$rhoMM ) )
  expect_identical( coef( flow_lin_sub2$eta ) , surf:::coef.svymstat( flow_lin_sub2$eta ) )
  expect_identical( coef( flow_lin_sub2$muij ) , surf:::coef.svymstat( flow_lin_sub2$muij ) )
  expect_identical( coef( flow_lin_sub2$pij ) , surf:::coef.svymstat( flow_lin_sub2$pij ) )
  expect_identical( coef( flow_rep_sub1$psi ) , surf:::coef.svymstat( flow_rep_sub1$psi ) )
  expect_identical( coef( flow_rep_sub1$rhoRR ) , surf:::coef.svymstat( flow_rep_sub1$rhoRR ) )
  expect_identical( coef( flow_rep_sub1$rhoMM ) , surf:::coef.svymstat( flow_rep_sub1$rhoMM ) )
  expect_identical( coef( flow_rep_sub1$eta ) , surf:::coef.svymstat( flow_rep_sub1$eta ) )
  expect_identical( coef( flow_rep_sub1$muij ) , surf:::coef.svymstat( flow_rep_sub1$muij ) )
  expect_identical( coef( flow_rep_sub1$pij ) , surf:::coef.svymstat( flow_rep_sub1$pij ) )
  expect_identical( coef( flow_rep_sub2$psi ) , surf:::coef.svymstat( flow_rep_sub2$psi ) )
  expect_identical( coef( flow_rep_sub2$rhoRR ) , surf:::coef.svymstat( flow_rep_sub2$rhoRR ) )
  expect_identical( coef( flow_rep_sub2$rhoMM ) , surf:::coef.svymstat( flow_rep_sub2$rhoMM ) )
  expect_identical( coef( flow_rep_sub2$eta ) , surf:::coef.svymstat( flow_rep_sub2$eta ) )
  expect_identical( coef( flow_rep_sub2$muij ) , surf:::coef.svymstat( flow_rep_sub2$muij ) )
  expect_identical( coef( flow_rep_sub2$pij ) , surf:::coef.svymstat( flow_rep_sub2$pij ) )

  # variances
  expect_identical( vcov( flow_lin_sub1$psi ) , survey:::vcov.svystat( flow_lin_sub1$psi ) )
  expect_identical( vcov( flow_lin_sub1$rhoRR ) , survey:::vcov.svystat( flow_lin_sub1$rhoRR ) )
  expect_identical( vcov( flow_lin_sub1$rhoMM ) , survey:::vcov.svystat( flow_lin_sub1$rhoMM ) )
  expect_identical( vcov( flow_lin_sub1$eta ) , survey:::vcov.svystat( flow_lin_sub1$eta ) )
  expect_identical( vcov( flow_lin_sub1$muij ) , surf:::vcov.svymstat( flow_lin_sub1$muij ) )
  expect_identical( vcov( flow_lin_sub1$pij ) , surf:::vcov.svymstat( flow_lin_sub1$pij ) )
  expect_identical( vcov( flow_lin_sub2$psi ) , survey:::vcov.svystat( flow_lin_sub2$psi ) )
  expect_identical( vcov( flow_lin_sub2$rhoRR ) , survey:::vcov.svystat( flow_lin_sub2$rhoRR ) )
  expect_identical( vcov( flow_lin_sub2$rhoMM ) , survey:::vcov.svystat( flow_lin_sub2$rhoMM ) )
  expect_identical( vcov( flow_lin_sub2$eta ) , survey:::vcov.svystat( flow_lin_sub2$eta ) )
  expect_identical( vcov( flow_lin_sub2$muij ) , surf:::vcov.svymstat( flow_lin_sub2$muij ) )
  expect_identical( vcov( flow_lin_sub2$pij ) , surf:::vcov.svymstat( flow_lin_sub2$pij ) )
  expect_identical( vcov( flow_rep_sub1$psi ) , survey:::vcov.svystat( flow_rep_sub1$psi ) )
  expect_identical( vcov( flow_rep_sub1$rhoRR ) , survey:::vcov.svystat( flow_rep_sub1$rhoRR ) )
  expect_identical( vcov( flow_rep_sub1$rhoMM ) , survey:::vcov.svystat( flow_rep_sub1$rhoMM ) )
  expect_identical( vcov( flow_rep_sub1$eta ) , survey:::vcov.svystat( flow_rep_sub1$eta ) )
  expect_identical( vcov( flow_rep_sub1$muij ) , surf:::vcov.svymstat( flow_rep_sub1$muij ) )
  expect_identical( vcov( flow_rep_sub1$pij ) , surf:::vcov.svymstat( flow_rep_sub1$pij ) )
  expect_identical( vcov( flow_rep_sub2$psi ) , survey:::vcov.svystat( flow_rep_sub2$psi ) )
  expect_identical( vcov( flow_rep_sub2$rhoRR ) , survey:::vcov.svystat( flow_rep_sub2$rhoRR ) )
  expect_identical( vcov( flow_rep_sub2$rhoMM ) , survey:::vcov.svystat( flow_rep_sub2$rhoMM ) )
  expect_identical( vcov( flow_rep_sub2$eta ) , survey:::vcov.svystat( flow_rep_sub2$eta ) )
  expect_identical( vcov( flow_rep_sub2$muij ) , surf:::vcov.svymstat( flow_rep_sub2$muij ) )
  expect_identical( vcov( flow_rep_sub2$pij ) , surf:::vcov.svymstat( flow_rep_sub2$pij ) )

  # standard errors
  expect_identical( SE( flow_lin_sub1$psi ) , survey:::SE.svystat( flow_lin_sub1$psi ) )
  expect_identical( SE( flow_lin_sub1$rhoRR ) , survey:::SE.svystat( flow_lin_sub1$rhoRR ) )
  expect_identical( SE( flow_lin_sub1$rhoMM ) , survey:::SE.svystat( flow_lin_sub1$rhoMM ) )
  expect_identical( SE( flow_lin_sub1$eta ) , survey:::SE.svystat( flow_lin_sub1$eta ) )
  expect_identical( SE( flow_lin_sub1$muij ) , surf:::SE.svymstat( flow_lin_sub1$muij ) )
  expect_identical( SE( flow_lin_sub1$pij ) , surf:::SE.svymstat( flow_lin_sub1$pij ) )
  expect_identical( SE( flow_lin_sub2$psi ) , survey:::SE.svystat( flow_lin_sub2$psi ) )
  expect_identical( SE( flow_lin_sub2$rhoRR ) , survey:::SE.svystat( flow_lin_sub2$rhoRR ) )
  expect_identical( SE( flow_lin_sub2$rhoMM ) , survey:::SE.svystat( flow_lin_sub2$rhoMM ) )
  expect_identical( SE( flow_lin_sub2$eta ) , survey:::SE.svystat( flow_lin_sub2$eta ) )
  expect_identical( SE( flow_lin_sub2$muij ) , surf:::SE.svymstat( flow_lin_sub2$muij ) )
  expect_identical( SE( flow_lin_sub2$pij ) , surf:::SE.svymstat( flow_lin_sub2$pij ) )
  expect_identical( SE( flow_rep_sub1$psi ) , survey:::SE.svystat( flow_rep_sub1$psi ) )
  expect_identical( SE( flow_rep_sub1$rhoRR ) , survey:::SE.svystat( flow_rep_sub1$rhoRR ) )
  expect_identical( SE( flow_rep_sub1$rhoMM ) , survey:::SE.svystat( flow_rep_sub1$rhoMM ) )
  expect_identical( SE( flow_rep_sub1$eta ) , survey:::SE.svystat( flow_rep_sub1$eta ) )
  expect_identical( SE( flow_rep_sub1$muij ) , surf:::SE.svymstat( flow_rep_sub1$muij ) )
  expect_identical( SE( flow_rep_sub1$pij ) , surf:::SE.svymstat( flow_rep_sub1$pij ) )
  expect_identical( SE( flow_rep_sub2$psi ) , survey:::SE.svystat( flow_rep_sub2$psi ) )
  expect_identical( SE( flow_rep_sub2$rhoRR ) , survey:::SE.svystat( flow_rep_sub2$rhoRR ) )
  expect_identical( SE( flow_rep_sub2$rhoMM ) , survey:::SE.svystat( flow_rep_sub2$rhoMM ) )
  expect_identical( SE( flow_rep_sub2$eta ) , survey:::SE.svystat( flow_rep_sub2$eta ) )
  expect_identical( SE( flow_rep_sub2$muij ) , surf:::SE.svymstat( flow_rep_sub2$muij ) )
  expect_identical( SE( flow_rep_sub2$pij ) , surf:::SE.svymstat( flow_rep_sub2$pij ) )

} )

# test against bias
test_that("compare point estimates vs population values",{

  expect_equivalent( coef( flow_lin_sub1$psi ) , psi_pop , tolerance = .1 )
  expect_equivalent( coef( flow_lin_sub2$psi ) , psi_pop , tolerance = .1 )
  expect_equivalent( coef( flow_lin_sub1$rhoRR ) , rhoRR_pop , tolerance = .1 )
  expect_equivalent( coef( flow_lin_sub2$rhoRR ) , rhoRR_pop , tolerance = .1 )
  expect_equivalent( coef( flow_lin_sub1$rhoMM ) , rhoMM_pop , tolerance = .1 )
  expect_equivalent( coef( flow_lin_sub2$rhoMM ) , rhoMM_pop , tolerance = .1 )
  expect_equivalent( coef( flow_lin_sub1$eta ) , eta_pop , tolerance = .1 )
  expect_equivalent( coef( flow_lin_sub2$eta ) , eta_pop , tolerance = .1 )
  expect_equivalent( coef( flow_lin_sub1$pij ) , pij_pop , tolerance = .1 )
  expect_equivalent( coef( flow_lin_sub2$muij ) , muij_pop / 2 , tolerance = .1 )

} )

# linearized vs replicate
test_that("compare linearized vs replicate: point estimates",{

  expect_identical( coef( flow_lin_sub1$psi ) , coef( flow_rep_sub1$psi ) )
  expect_identical( coef( flow_lin_sub1$rhoRR ) , coef( flow_rep_sub1$rhoRR ) )
  expect_identical( coef( flow_lin_sub1$rhoMM ) , coef( flow_rep_sub1$rhoMM ) )
  expect_identical( coef( flow_lin_sub1$eta ) , coef( flow_rep_sub1$eta ) )
  expect_identical( coef( flow_lin_sub1$pij ) , coef( flow_rep_sub1$pij ) )
  expect_identical( coef( flow_lin_sub1$muij ) , coef( flow_rep_sub1$muij ) )

  expect_identical( coef( flow_lin_sub2$psi ) , coef( flow_rep_sub2$psi ) )
  expect_identical( coef( flow_lin_sub2$rhoRR ) , coef( flow_rep_sub2$rhoRR ) )
  expect_identical( coef( flow_lin_sub2$rhoMM ) , coef( flow_rep_sub2$rhoMM ) )
  expect_identical( coef( flow_lin_sub2$eta ) , coef( flow_rep_sub2$eta ) )
  expect_identical( coef( flow_lin_sub2$pij ) , coef( flow_rep_sub2$pij ) )
  expect_identical( coef( flow_lin_sub2$muij ) , coef( flow_rep_sub2$muij ) )

} )
test_that("compare linearized vs replicate: standard errors",{

  expect_equivalent( SE( flow_lin_sub1$psi ) , SE( flow_rep_sub1$psi ) , tolerance = .1 )
  expect_equivalent( SE( flow_lin_sub1$rhoRR ) , SE( flow_rep_sub1$rhoRR ) , tolerance = .1 )
  expect_equivalent( SE( flow_lin_sub1$rhoMM ) , SE( flow_rep_sub1$rhoMM ) , tolerance = .1 )
  expect_equivalent( SE( flow_lin_sub1$eta ) , SE( flow_rep_sub1$eta ) , tolerance = .1 )
  expect_equivalent( SE( flow_lin_sub1$pij ) , SE( flow_rep_sub1$pij ) , tolerance = .1 )
  expect_equivalent( SE( flow_lin_sub1$muij ) , SE( flow_rep_sub1$muij ) , tolerance = .4 )

  expect_equivalent( SE( flow_lin_sub2$psi ) , SE( flow_rep_sub2$psi ) , tolerance = .1 )
  expect_equivalent( SE( flow_lin_sub2$rhoRR ) , SE( flow_rep_sub2$rhoRR ) , tolerance = .1 )
  expect_equivalent( SE( flow_lin_sub2$rhoMM ) , SE( flow_rep_sub2$rhoMM ) , tolerance = .1 )
  expect_equivalent( SE( flow_lin_sub2$eta ) , SE( flow_rep_sub2$eta ) , tolerance = .1 )
  expect_equivalent( SE( flow_lin_sub2$pij ) , SE( flow_rep_sub2$pij ) , tolerance = .1 )
  expect_equivalent( SE( flow_lin_sub2$muij ) , SE( flow_rep_sub2$muij ) , tolerance = .4 )

} )
