context("svyflow estimates with nonresponse")

# set random seed
set.seed(123)

# create pseudo-population
N = 1e5
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

# adds non-response
pop_frame[ as.logical( rbinom(N,1,.2) ) , 1 ] <- NA
pop_frame[ !is.na( pop_frame[,1] ) & as.logical( rbinom(N,1,.1) ) , 2 ] <- NA
pop_frame[ is.na( pop_frame[,2] ) & as.logical( rbinom(N,1,.7) ) , 1 ] <- NA

# extract sample
smp_frame <- pop_frame[ as.logical( sampling::srswor( n , N ) ) , ]
rm( pop_frame )

# splits data
df0 <- smp_frame[ , c(1,3) ]
colnames(df0) <- c("v0","v1")
df0$prob <- n/N
df1 <- smp_frame[ , c(2,4) ]
colnames(df1) <- c("v0","v1")

# load libraries
library( surf )
library( testthat )

# build designs
flowdes_srs <-
  sfydesign( ids = ~0 ,
             fpc = ~ 1/prob ,
             data = list( df0 , df1 ) ,
             nest = TRUE )

# create replicate designs
flowdes_srs_rep <- as.surfrdesign( flowdes_srs , type = "bootstrap" , replicate = 50 )

# estimate flows
gflow_srs_lin <- svyflow( ~v0 , flowdes_srs , na.rm = TRUE )
nflow_srs_lin <- svyflow( ~v0 , flowdes_srs , flow.type = "net" , na.rm = TRUE )
gflow_srs_rep <- svyflow( ~v0 , flowdes_srs_rep , na.rm = TRUE )
nflow_srs_rep <- svyflow( ~v0 , flowdes_srs_rep , flow.type = "net" , na.rm = TRUE )
gflow_srs_lin_na <- svyflow( ~v0 , flowdes_srs )
nflow_srs_lin_na <- svyflow( ~v0 , flowdes_srs , flow.type = "net" )
gflow_srs_rep_na <- svyflow( ~v0 , flowdes_srs_rep )
nflow_srs_rep_na <- svyflow( ~v0 , flowdes_srs_rep , flow.type = "net" )

# test extraction of associated measures
test_that( "extraction of estimates" , {

  # point estimates
  expect_identical( coef( gflow_srs_lin ) , surf:::coef.flowstat( gflow_srs_lin ) )
  expect_identical( coef( gflow_srs_rep ) , surf:::coef.flowstat( gflow_srs_rep ) )
  expect_identical( coef( nflow_srs_lin ) , surf:::coef.flowstat( nflow_srs_lin ) )
  expect_identical( coef( nflow_srs_rep ) , surf:::coef.flowstat( nflow_srs_rep ) )

  # variances
  expect_identical( vcov( gflow_srs_lin ) , surf:::vcov.flowstat( gflow_srs_lin ) )
  expect_identical( vcov( gflow_srs_rep ) , surf:::vcov.flowstat( gflow_srs_rep ) )
  expect_identical( vcov( nflow_srs_lin ) , surf:::vcov.flowstat( nflow_srs_lin ) )
  expect_identical( vcov( nflow_srs_rep ) , surf:::vcov.flowstat( nflow_srs_rep ) )

  # standard errors
  expect_identical( SE( gflow_srs_lin ) , surf:::SE.flowstat( gflow_srs_lin ) )
  expect_identical( SE( gflow_srs_rep ) , surf:::SE.flowstat( gflow_srs_rep ) )
  expect_identical( SE( nflow_srs_lin ) , surf:::SE.flowstat( nflow_srs_lin ) )
  expect_identical( SE( nflow_srs_rep ) , surf:::SE.flowstat( nflow_srs_rep ) )

  # coefficients of variation
  expect_identical( cv( gflow_srs_lin ) , surf:::cv.flowstat( gflow_srs_lin ) )
  expect_identical( cv( gflow_srs_rep ) , surf:::cv.flowstat( gflow_srs_rep ) )
  expect_identical( cv( nflow_srs_lin ) , surf:::cv.flowstat( nflow_srs_lin ) )
  expect_identical( cv( nflow_srs_rep ) , surf:::cv.flowstat( nflow_srs_rep ) )

} )

# test against bias
test_that("compare point estimates vs population values",{
  expect_equivalent( coef( gflow_srs_lin ) , muij_pop , tolerance = .05 )
  expect_equivalent( coef( gflow_srs_rep ) , muij_pop , tolerance = .05 )
  expect_equivalent( coef( nflow_srs_lin ) , nipij_pop , tolerance = .05 )
  expect_equivalent( coef( nflow_srs_rep ) , nipij_pop , tolerance = .05 )
} )

# linearized vs replicate
test_that("compare linearized vs replicate: point estimates",{
  expect_identical( coef( gflow_srs_lin ) , coef( gflow_srs_rep ) )
  expect_identical( coef( nflow_srs_lin ) , coef( nflow_srs_rep ) )
} )
test_that("compare linearized vs replicate: standard errors",{
  expect_equal( SE( gflow_srs_lin ) , SE( gflow_srs_rep ) , tolerance = .4 )
  expect_equal( SE( nflow_srs_lin ) , SE( nflow_srs_rep ) , tolerance = .4 )
} )
test_that("compare linearized vs replicate: coefficients of variation",{
  expect_equal( cv( gflow_srs_lin ) , cv( gflow_srs_rep ) , tolerance = .4 )
  expect_equal( cv( nflow_srs_lin ) , cv( nflow_srs_rep ) , tolerance = .4 )
} )

# # check for consistency across versions
# test_that( "version-consistency tests" , {
#
#   # point estimates
#   verify_output( "output/gflow_srs_lin_coef.txt" , coef( gflow_srs_lin ) )
#   verify_output( "output/gflow_srs_rep_coef.txt" , coef( gflow_srs_rep ) )
#   verify_output( "output/nflow_srs_lin_coef.txt" , coef( nflow_srs_lin ) )
#   verify_output( "output/nflow_srs_rep_coef.txt" , coef( nflow_srs_rep ) )
#
#   # variances
#   verify_output( "output/gflow_srs_lin_vcov.txt" , vcov( gflow_srs_lin ) )
#   # verify_output( "output/gflow_srs_rep_vcov.txt" , vcov( gflow_srs_rep ) )
#   verify_output( "output/nflow_srs_lin_vcov.txt" , vcov( nflow_srs_lin ) )
#   # verify_output( "output/nflow_srs_rep_vcov.txt" , vcov( nflow_srs_rep ) )
#
#   # standard errors
#   verify_output( "output/gflow_srs_lin_SE.txt" , SE( gflow_srs_lin ) )
#   # verify_output( "output/gflow_srs_rep_SE.txt" , SE( gflow_srs_rep ) )
#   verify_output( "output/nflow_srs_lin_SE.txt" , SE( nflow_srs_lin ) )
#   # verify_output( "output/nflow_srs_rep_SE.txt" , SE( nflow_srs_rep ) )
#
#   # coefficients of variation
#   verify_output( "output/gflow_srs_lin_cv.txt" , cv( gflow_srs_lin ) )
#   # verify_output( "output/gflow_srs_rep_cv.txt" , cv( gflow_srs_rep ) )
#   verify_output( "output/nflow_srs_lin_cv.txt" , cv( nflow_srs_lin ) )
#   # verify_output( "output/nflow_srs_rep_cv.txt" , cv( nflow_srs_rep ) )
#
# } )
