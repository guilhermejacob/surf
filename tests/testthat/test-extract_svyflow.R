context("flow output")

# build designs
flowdes <-
  sfydesign( ids = ~0 ,
             probs = ~ prob ,
             data = list( dfa0 , dfa1 ) ,
             nest = TRUE )

# create replicate design
flowdes_rep <- as.surfrdesign( flowdes , type = "bootstrap" , replicate = 50 )

# estimate flows
gflow_lin <- svyflow( ~v0 , flowdes )
nflow_lin <- svyflow( ~v0 , flowdes , flow.type = "net" )
gflow_rep <- svyflow( ~v0 , flowdes_rep )
nflow_rep <- svyflow( ~v0 , flowdes_rep , flow.type = "net" )

test_that("outputs",{
  expect_is(gflow_lin,"flowstat")
  expect_is(gflow_rep,"flowstat")
  expect_is(nflow_rep,"flowstat")
  expect_is(nflow_lin,"flowstat")
  expect_is(coef(gflow_lin),"matrix")
  expect_is(coef(gflow_rep),"matrix")
  expect_is(coef(nflow_rep),"matrix")
  expect_is(coef(nflow_lin),"matrix")
  expect_equal(sum(coef(nflow_lin)),1)
  expect_equal(sum(coef(nflow_rep)),1)
  expect_identical(coef(gflow_lin), coef(gflow_rep))
  expect_identical(coef(nflow_lin), coef(nflow_lin))
  expect_is(SE(gflow_lin),"table")
  expect_is(SE(nflow_lin),"table")
  expect_is(SE(gflow_rep),"table")
  expect_is(SE(nflow_rep),"table")
})

