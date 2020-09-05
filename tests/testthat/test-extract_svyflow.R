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

# test outputs
test_that("outputs",{

  expect_is( flow_srs_lin , "flowstat" )
  # expect_is( flow_srs_rep , "flowstat" )

  expect_is( flow_srs_lin$psi , "svystat" )
  expect_is( flow_srs_lin$rho , "svystat" )
  expect_is( flow_srs_lin$ta , "svystat" )
  expect_is( flow_srs_lin$eta , "svystat" )
  expect_is( flow_srs_lin$pij , "svymstat" )
  expect_is( flow_srs_lin$muij , "svymstat" )
  # expect_is( flow_srs_rep$psi , "svystat" )
  # expect_is( flow_srs_rep$tau , "svystat" )
  # expect_is( flow_srs_rep$tau , "svystat" )
  # expect_is( flow_srs_rep$eta , "svystat" )
  # expect_is( flow_srs_rep$pij , "svymstat" )
  # expect_is( flow_srs_rep$muij , "svymstat" )

  expect_equivalent( sum( coef( flow_srs_lin$eta ) ) , 1 )
  expect_equivalent( rowSums( coef( flow_srs_lin$pij ) ) , c(1,1,1) )
  # expect_equivalent( sum( coef( flow_srs_rep$eta ) ) , 1 )
  # expect_equivalent( rowSums( coef( flow_srs_rep$pij ) ) , c(1,1,1) )

} )
