context("flow output")

# random number generator seed
set.seed( 123 )

# population size
N = as.integer( 10^5 )

# sample size
n = as.integer(10^4)

# superopopulation model parameters
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

# cria objeto de desenho amostral
df0 <- smp_df[ , -2 , drop = FALSE ]
df1 <- smp_df[ , 2 , drop = FALSE ]
colnames( df0 )[1] <- "y"
colnames( df1 )[1] <- "y"

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
flow_srs_rep <- svyflow( ~y , flowdes_srs_rep )

# test outputs
test_that("outputs",{

  expect_is( flow_srs_lin , "flowstat" )
  expect_is( flow_srs_rep , "flowstat" )

  expect_is( flow_srs_lin$psi , "svystat" )
  expect_is( flow_srs_lin$rhoRR , "svystat" )
  expect_is( flow_srs_lin$rhoMM , "svystat" )
  expect_is( flow_srs_lin$eta , "svystat" )
  expect_is( flow_srs_lin$pij , "svymstat" )
  expect_is( flow_srs_lin$muij , "svymstat" )
  expect_is( flow_srs_rep$psi , "svystat" )
  expect_is( flow_srs_rep$rhoRR , "svystat" )
  expect_is( flow_srs_rep$rhoMM , "svystat" )
  expect_is( flow_srs_rep$eta , "svystat" )
  expect_is( flow_srs_rep$pij , "svymstat" )
  expect_is( flow_srs_rep$muij , "svymstat" )

  expect_equivalent( sum( coef( flow_srs_lin$eta ) ) , 1 )
  expect_equivalent( rowSums( coef( flow_srs_lin$pij ) ) , c(1,1,1) )
  expect_equivalent( sum( coef( flow_srs_rep$eta ) ) , 1 )
  expect_equivalent( rowSums( coef( flow_srs_rep$pij ) ) , c(1,1,1) )

} )
