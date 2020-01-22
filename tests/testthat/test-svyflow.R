context("flow estimates")

# define seed for random numbers
set.seed(123)

# build designs
flowdes_srs <-
  sfydesign( ids = ~0 ,
             fpc = ~ 1/prob ,
             data = list( dfa0 , dfa1 ) ,
             nest = TRUE )
flowdes_strat <-
  sfydesign( ids = ~0 ,
             probs = ~ prob ,
             strata = ~ stratum ,
             data = list( surf:::dfa0_strat , surf:::dfa1_strat ) ,
             nest = TRUE )

# create replicate designs
flowdes_srs_rep <- as.surfrdesign( flowdes_srs , type = "bootstrap" , replicate = 50 )
flowdes_strat_rep <- as.surfrdesign( flowdes_strat , type = "bootstrap" , replicate = 50 )

# estimate flows
gflow_srs_lin <- svyflow( ~v0 , flowdes_srs )
gflow_strat_lin <- svyflow( ~v0 , flowdes_strat )
nflow_srs_lin <- svyflow( ~v0 , flowdes_srs , flow.type = "net" )
nflow_strat_lin <- svyflow( ~v0 , flowdes_strat , flow.type = "net" )
gflow_srs_rep <- svyflow( ~v0 , flowdes_srs_rep )
gflow_strat_rep <- svyflow( ~v0 , flowdes_strat_rep )
nflow_srs_rep <- svyflow( ~v0 , flowdes_srs_rep , flow.type = "net" )
nflow_strat_rep <- svyflow( ~v0 , flowdes_strat_rep , flow.type = "net" )

# test extraction of associated measures
test_that( "extraction of estimates" , {

  # point estimates
  expect_identical( coef( gflow_srs_lin ) , surf:::coef.flowstat( gflow_srs_lin ) )
  expect_identical( coef( gflow_strat_lin ) , surf:::coef.flowstat( gflow_strat_lin ) )
  expect_identical( coef( gflow_srs_rep ) , surf:::coef.flowstat( gflow_srs_rep ) )
  expect_identical( coef( gflow_strat_rep ) , surf:::coef.flowstat( gflow_strat_rep ) )
  expect_identical( coef( nflow_srs_lin ) , surf:::coef.flowstat( nflow_srs_lin ) )
  expect_identical( coef( nflow_strat_lin ) , surf:::coef.flowstat( nflow_strat_lin ) )
  expect_identical( coef( nflow_srs_rep ) , surf:::coef.flowstat( nflow_srs_rep ) )
  expect_identical( coef( nflow_strat_rep ) , surf:::coef.flowstat( nflow_strat_rep ) )

  # variances
  expect_identical( vcov( gflow_srs_lin ) , surf:::vcov.flowstat( gflow_srs_lin ) )
  expect_identical( vcov( gflow_strat_lin ) , surf:::vcov.flowstat( gflow_strat_lin ) )
  expect_identical( vcov( gflow_srs_rep ) , surf:::vcov.flowstat( gflow_srs_rep ) )
  expect_identical( vcov( gflow_strat_rep ) , surf:::vcov.flowstat( gflow_strat_rep ) )
  expect_identical( vcov( nflow_srs_lin ) , surf:::vcov.flowstat( nflow_srs_lin ) )
  expect_identical( vcov( nflow_strat_lin ) , surf:::vcov.flowstat( nflow_strat_lin ) )
  expect_identical( vcov( nflow_srs_rep ) , surf:::vcov.flowstat( nflow_srs_rep ) )
  expect_identical( vcov( nflow_strat_rep ) , surf:::vcov.flowstat( nflow_strat_rep ) )

  # standard errors
  expect_identical( SE( gflow_srs_lin ) , surf:::SE.flowstat( gflow_srs_lin ) )
  expect_identical( SE( gflow_strat_lin ) , surf:::SE.flowstat( gflow_strat_lin ) )
  expect_identical( SE( gflow_srs_rep ) , surf:::SE.flowstat( gflow_srs_rep ) )
  expect_identical( SE( gflow_strat_rep ) , surf:::SE.flowstat( gflow_strat_rep ) )
  expect_identical( SE( nflow_srs_lin ) , surf:::SE.flowstat( nflow_srs_lin ) )
  expect_identical( SE( nflow_strat_lin ) , surf:::SE.flowstat( nflow_strat_lin ) )
  expect_identical( SE( nflow_srs_rep ) , surf:::SE.flowstat( nflow_srs_rep ) )
  expect_identical( SE( nflow_strat_rep ) , surf:::SE.flowstat( nflow_strat_rep ) )

  # coefficients of variation
  expect_identical( cv( gflow_srs_lin ) , surf:::cv.flowstat( gflow_srs_lin ) )
  expect_identical( cv( gflow_strat_lin ) , surf:::cv.flowstat( gflow_strat_lin ) )
  expect_identical( cv( gflow_srs_rep ) , surf:::cv.flowstat( gflow_srs_rep ) )
  expect_identical( cv( gflow_strat_rep ) , surf:::cv.flowstat( gflow_strat_rep ) )
  expect_identical( cv( nflow_srs_lin ) , surf:::cv.flowstat( nflow_srs_lin ) )
  expect_identical( cv( nflow_strat_lin ) , surf:::cv.flowstat( nflow_strat_lin ) )
  expect_identical( cv( nflow_srs_rep ) , surf:::cv.flowstat( nflow_srs_rep ) )
  expect_identical( cv( nflow_strat_rep ) , surf:::cv.flowstat( nflow_strat_rep ) )

} )

# check for consistency across versions
test_that( "version-consistency tests" , {

  # point estimates
  verify_output( "gflow_srs_lin_coef.txt" , coef( gflow_srs_lin ) )
  verify_output( "gflow_strat_lin_coef.txt" , coef( gflow_strat_lin ) )
  verify_output( "gflow_srs_rep_coef.txt" , coef( gflow_srs_rep ) )
  verify_output( "gflow_strat_rep_coef.txt" , coef( gflow_strat_rep ) )
  verify_output( "nflow_srs_lin_coef.txt" , coef( nflow_srs_lin ) )
  verify_output( "nflow_strat_lin_coef.txt" , coef( nflow_strat_lin ) )
  verify_output( "nflow_srs_rep_coef.txt" , coef( nflow_srs_rep ) )
  verify_output( "nflow_strat_rep_coef.txt" , coef( nflow_strat_rep ) )

  # variances
  verify_output( "gflow_srs_lin_vcov.txt" , vcov( gflow_srs_lin ) )
  verify_output( "gflow_strat_lin_vcov.txt" , vcov( gflow_strat_lin ) )
  verify_output( "gflow_srs_rep_vcov.txt" , vcov( gflow_srs_rep ) )
  verify_output( "gflow_strat_rep_vcov.txt" , vcov( gflow_strat_rep ) )
  verify_output( "nflow_srs_lin_vcov.txt" , vcov( nflow_srs_lin ) )
  verify_output( "nflow_strat_lin_vcov.txt" , vcov( nflow_strat_lin ) )
  verify_output( "nflow_srs_rep_vcov.txt" , vcov( nflow_srs_rep ) )
  verify_output( "nflow_strat_rep_vcov.txt" , vcov( nflow_strat_rep ) )

  # standard errors
  verify_output( "gflow_srs_lin_SE.txt" , SE( gflow_srs_lin ) )
  verify_output( "gflow_strat_lin_SE.txt" , SE( gflow_strat_lin ) )
  verify_output( "gflow_srs_rep_SE.txt" , SE( gflow_srs_rep ) )
  verify_output( "gflow_strat_rep_SE.txt" , SE( gflow_strat_rep ) )
  verify_output( "nflow_srs_lin_SE.txt" , SE( nflow_srs_lin ) )
  verify_output( "nflow_strat_lin_SE.txt" , SE( nflow_strat_lin ) )
  verify_output( "nflow_srs_rep_SE.txt" , SE( nflow_srs_rep ) )
  verify_output( "nflow_strat_rep_SE.txt" , SE( nflow_strat_rep ) )

  # coefficients of variation
  verify_output( "gflow_srs_lin_cv.txt" , cv( gflow_srs_lin ) )
  verify_output( "gflow_strat_lin_cv.txt" , cv( gflow_strat_lin ) )
  verify_output( "gflow_srs_rep_cv.txt" , cv( gflow_srs_rep ) )
  verify_output( "gflow_strat_rep_cv.txt" , cv( gflow_strat_rep ) )
  verify_output( "nflow_srs_lin_cv.txt" , cv( nflow_srs_lin ) )
  verify_output( "nflow_strat_lin_cv.txt" , cv( nflow_strat_lin ) )
  verify_output( "nflow_srs_rep_cv.txt" , cv( nflow_srs_rep ) )
  verify_output( "nflow_strat_rep_cv.txt" , cv( nflow_strat_rep ) )

} )

# test against bias
test_that("compare point estimates vs population values",{
  expect_equivalent( coef( gflow_srs_lin ) , surf:::muij_pop , tolerance = .05 )
  expect_equivalent( coef( gflow_srs_rep ) , surf:::muij_pop , tolerance = .05 )
  expect_equivalent( coef( gflow_strat_lin ) , surf:::muij_pop , tolerance = .05 )
  expect_equivalent( coef( gflow_strat_rep ) , surf:::muij_pop , tolerance = .05 )
  expect_equivalent( coef( nflow_srs_lin ) , surf:::nipij_pop , tolerance = .05 )
  expect_equivalent( coef( nflow_srs_rep ) , surf:::nipij_pop , tolerance = .05 )
  expect_equivalent( coef( nflow_strat_lin ) , surf:::nipij_pop , tolerance = .05 )
  expect_equivalent( coef( nflow_strat_rep ) , surf:::nipij_pop , tolerance = .05 )
} )

# linearized vs replicate
test_that("compare linearized vs replicate: point estimates",{
  expect_identical( coef( gflow_srs_lin ) , coef( gflow_srs_rep ) )
  expect_identical( coef( gflow_strat_lin ) , coef( gflow_strat_rep ) )
  expect_identical( coef( nflow_srs_lin ) , coef( nflow_srs_rep ) )
  expect_identical( coef( nflow_strat_lin ) , coef( nflow_strat_rep ) )
} )
test_that("compare linearized vs replicate: standard errors",{
  expect_equal( SE( gflow_srs_lin ) , SE( gflow_srs_rep ) , tolerance = .4 )
  expect_equal( SE( gflow_strat_lin ) , SE( gflow_strat_rep ) , tolerance = .4 )
  expect_equal( SE( nflow_srs_lin ) , SE( nflow_srs_rep ) , tolerance = .4 )
  expect_equal( SE( nflow_strat_lin ) , SE( nflow_strat_rep ) , tolerance = .4 )
} )
test_that("compare linearized vs replicate: coefficients of variation",{
  expect_equal( cv( gflow_srs_lin ) , cv( gflow_srs_rep ) , tolerance = .3 )
  expect_equal( cv( gflow_strat_lin ) , cv( gflow_strat_rep ) , tolerance = .3 )
  expect_equal( cv( nflow_srs_lin ) , cv( nflow_srs_rep ) , tolerance = .3 )
  expect_equal( cv( nflow_strat_lin ) , cv( nflow_strat_rep ) , tolerance = .3 )
} )
