context("updates on survflow.design")

# build designs
flowdes <-
  sfydesign( ids = ~0 ,
             fpc = ~ 1/prob ,
             data = list( dfa0 , dfa1 ) ,
             nest = TRUE )

# create replicate designs
flowdes_rep <- as.surfrdesign( flowdes , type = "bootstrap" , replicate = 50 )

# update subsets
flowdes_trans <- update( flowdes , idstatus = factor( as.numeric(v0) , labels = LETTERS[1:3] , levels = 1:3 ) )
flowdes_rep_trans <- update( flowdes_rep , idstatus = factor( as.numeric(v0) , labels = LETTERS[1:3] , levels = 1:3 ) )

# estimate flows
gflow_lin <- svyflow( ~v0 , flowdes )
gflow_lin_trans <- svyflow( ~idstatus , flowdes_trans )
gflow_rep <- svyflow( ~v0 , flowdes_rep )
gflow_rep_trans <- svyflow( ~idstatus , flowdes_rep_trans )
nflow_lin <- svyflow( ~v0 , flowdes , flow.type = "net" )
nflow_lin_trans <- svyflow( ~idstatus , flowdes_trans , flow.type = "net" )
nflow_rep <- svyflow( ~v0 , flowdes_rep , flow.type = "net" )
nflow_rep_trans <- svyflow( ~idstatus , flowdes_rep_trans , flow.type = "net" )

# linearized vs replicate
test_that("compare linearized vs replicate: estimates",{
  expect_identical( coef( gflow_lin_trans ) , coef( gflow_rep_trans ) )
  expect_equal( SE( gflow_lin_trans ) , SE( gflow_rep_trans ) , tolerance = .3 )
  expect_equal( cv( gflow_lin_trans ) , cv( gflow_rep_trans ) , tolerance = .3 )
  expect_identical( coef( nflow_lin_trans ) , coef( nflow_rep_trans ) )
  expect_equal( SE( nflow_lin_trans ) , SE( nflow_rep_trans ) , tolerance = .3 )
  expect_equal( cv( nflow_lin_trans ) , cv( nflow_rep_trans ) , tolerance = .3 )
} )

# untransformed vs transformed
test_that("compare transformed vs unstransformed",{

  expect_failure( expect_identical( coef( gflow_lin ) , coef( gflow_lin_trans ) ) )
  expect_identical( matrix( coef( gflow_lin ) ) , matrix( coef( gflow_lin_trans ) ) )
  expect_identical( matrix( SE( gflow_lin ) ) , matrix( SE( gflow_lin_trans ) ) )
  expect_identical( matrix( cv( gflow_lin ) ) , matrix( cv( gflow_lin_trans ) ) )
  expect_failure( expect_identical( coef( nflow_lin ) , coef( nflow_lin_trans ) ) )
  expect_identical( matrix( coef( nflow_lin ) ) , matrix( coef( nflow_lin_trans ) ) )
  expect_identical( matrix( SE( nflow_lin ) ) , matrix( SE( nflow_lin_trans ) ) )
  expect_identical( matrix( cv( nflow_lin ) ) , matrix( cv( nflow_lin_trans ) ) )

} )
