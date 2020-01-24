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
gflow_lin <- svyflow( ~v0 , flowdes , na.rm = TRUE )
nflow_lin <- svyflow( ~v0 , flowdes , flow.type = "net" , na.rm = TRUE )
gflow_rep <- svyflow( ~v0 , flowdes_rep , na.rm = TRUE )
nflow_rep <- svyflow( ~v0 , flowdes_rep , flow.type = "net" , na.rm = TRUE )
gflow_lin_na <- svyflow( ~v0 , flowdes )
nflow_lin_na <- svyflow( ~v0 , flowdes , flow.type = "net" )
gflow_rep_na <- svyflow( ~v0 , flowdes_rep )
nflow_rep_na <- svyflow( ~v0 , flowdes_rep , flow.type = "net" )

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
  expect_is(SE(gflow_lin),"matrix")
  expect_is(SE(nflow_lin),"matrix")
  expect_is(SE(gflow_rep),"matrix")
  expect_is(SE(nflow_rep),"matrix")
})

test_that("outputs with NA",{
  expect_is(gflow_lin_na,"flowstat")
  expect_is(gflow_rep_na,"flowstat")
  expect_is(nflow_rep_na,"flowstat")
  expect_is(nflow_lin_na,"flowstat")
  expect_is(coef(gflow_lin_na),"matrix")
  expect_is(coef(gflow_rep_na),"matrix")
  expect_is(coef(nflow_rep_na),"matrix")
  expect_is(coef(nflow_lin_na),"matrix")
  expect_equal(sum(coef(nflow_lin_na)),as.numeric(NA))
  expect_equal(sum(coef(nflow_rep_na)),as.numeric(NA))
  expect_equal(sum(SE(nflow_lin_na)),as.numeric(NA))
  expect_equal(sum(SE(nflow_rep_na)),as.numeric(NA))
  expect_identical(coef(gflow_lin_na), coef(gflow_rep_na))
  expect_identical(coef(nflow_lin_na), coef(nflow_lin_na))
  expect_is(SE(gflow_lin_na),"matrix")
  expect_is(SE(nflow_lin_na),"matrix")
  expect_is(SE(gflow_rep_na),"matrix")
  expect_is(SE(nflow_rep_na),"matrix")
})
