context("svyflow with full response")

# set random seed
set.seed(123)

# create pseudo-population
N = 10^5
n = round(N *.05)
eta_pop <- c(.9,.05,.05)
pij_pop <- matrix( c(0.8, 0.15, 0.05, 0.3, 0.6, 0.1, 0.1, 0.1, 0.8) , ncol = 3 , nrow = 3 , byrow = T )
nipij_pop <- sweep( pij_pop , 1 , eta_pop , "*" )
muij_pop <- nipij_pop * N
pop <- as.numeric( muij_pop )
pop <- data.frame( k_ij = do.call( c , sapply( seq_along( pop ), function(z) rep( z , pop[z] ) ) ) )
pop_frame <- expand.grid( data.frame( v0 = seq_len( nrow(nipij_pop) ) , v1 = seq_len( nrow(nipij_pop) ) ) )
pop_frame <- cbind( k_ij = seq_len( nrow( pop_frame ) ) , pop_frame )
pop_frame <- merge( pop , pop_frame )[ , -1 ] ; rm( pop )

# create covariates
pop_frame$v2 <- ifelse( pop_frame$v0 == 1 , rbinom( nrow( pop_frame ) , 1 , .7 ) , rbinom( nrow( pop_frame ) , 1 , .3 ) ) + 1
pop_frame$v3 <- ifelse( pop_frame$v0 == 1 , rbinom( nrow( pop_frame ) , 1 , .7 ) , rbinom( nrow( pop_frame ) , 1 , .3 ) ) + 1
pop_frame <- pop_frame[ order( pop_frame$v1 ) , ]
pop_frame <- data.frame( apply(pop_frame[,] , 2 , factor) , stringsAsFactors = TRUE )

# extract sample
smp_frame <- pop_frame[ as.logical( sampling::srswor( n , N ) ) , ]
# rm( pop_frame )

# splits data
dfa0 <- smp_frame[ , c(1,3) ]
colnames(dfa0) <- c("v0","v1")
dfa0$prob <- n/N
dfa1 <- smp_frame[ , c(2,4) ]
colnames(dfa1) <- c("v0","v1")

# build designs
flowdes <-
  sfydesign( ids = ~0 ,
             fpc = ~ 1/prob ,
             data = list( dfa0 , dfa1 ) ,
             nest = TRUE )

# create replicate designs
flowdes_rep <- as.surfrdesign( flowdes , type = "bootstrap" , replicate = 50 )

# estimate flows
gflow_lin <- svyflow( ~v0 , flowdes )
nflow_lin <- svyflow( ~v0 , flowdes , flow.type = "net" )
gflow_rep <- svyflow( ~v0 , flowdes_rep )
nflow_rep <- svyflow( ~v0 , flowdes_rep , flow.type = "net" )

# test extraction of associated measures
test_that( "extraction of estimates" , {

  # point estimates
  expect_identical( coef( gflow_lin ) , surf:::coef.flowstat( gflow_lin ) )
  expect_identical( coef( gflow_rep ) , surf:::coef.flowstat( gflow_rep ) )
  expect_identical( coef( nflow_lin ) , surf:::coef.flowstat( nflow_lin ) )
  expect_identical( coef( nflow_rep ) , surf:::coef.flowstat( nflow_rep ) )

  # variances
  expect_identical( vcov( gflow_lin ) , surf:::vcov.flowstat( gflow_lin ) )
  expect_identical( vcov( gflow_rep ) , surf:::vcov.flowstat( gflow_rep ) )
  expect_identical( vcov( nflow_lin ) , surf:::vcov.flowstat( nflow_lin ) )
  expect_identical( vcov( nflow_rep ) , surf:::vcov.flowstat( nflow_rep ) )

  # standard errors
  expect_identical( SE( gflow_lin ) , surf:::SE.flowstat( gflow_lin ) )
  expect_identical( SE( gflow_rep ) , surf:::SE.flowstat( gflow_rep ) )
  expect_identical( SE( nflow_lin ) , surf:::SE.flowstat( nflow_lin ) )
  expect_identical( SE( nflow_rep ) , surf:::SE.flowstat( nflow_rep ) )

  # coefficients of variation
  expect_identical( cv( gflow_lin ) , surf:::cv.flowstat( gflow_lin ) )
  expect_identical( cv( gflow_rep ) , surf:::cv.flowstat( gflow_rep ) )
  expect_identical( cv( nflow_lin ) , surf:::cv.flowstat( nflow_lin ) )
  expect_identical( cv( nflow_rep ) , surf:::cv.flowstat( nflow_rep ) )

} )

# test against bias
test_that("compare point estimates vs population values",{
  expect_equivalent( coef( gflow_lin ) , surf:::muij_pop , tolerance = .05 )
  expect_equivalent( coef( gflow_rep ) , surf:::muij_pop , tolerance = .05 )
  expect_equivalent( coef( nflow_lin ) , surf:::nipij_pop , tolerance = .05 )
  expect_equivalent( coef( nflow_rep ) , surf:::nipij_pop , tolerance = .05 )
} )

# linearized vs replicate
test_that("compare linearized vs replicate: point estimates",{
  expect_identical( coef( gflow_lin ) , coef( gflow_rep ) )
  expect_identical( coef( nflow_lin ) , coef( nflow_rep ) )
} )
test_that("compare linearized vs replicate: standard errors",{
  expect_equal( SE( gflow_lin ) , SE( gflow_rep ) , tolerance = .4 )
  expect_equal( SE( nflow_lin ) , SE( nflow_rep ) , tolerance = .4 )
} )
test_that("compare linearized vs replicate: coefficients of variation",{
  expect_equal( cv( gflow_lin ) , cv( gflow_rep ) , tolerance = .3 )
  expect_equal( cv( nflow_lin ) , cv( nflow_rep ) , tolerance = .3 )
} )
