context("subsets on survflow.design")

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

# splits data
popdata <- as.matrix( pop_frame[ , c(1,2) ] )
colnames(popdata) <- paste0( "round" , seq_len( ncol(popdata) ) - 1 , ":v0" )

# cross tabulation
muij_pop1 <- xtabs( ~. , popdata , subset = pop_frame$v3 == 1 )
muij_pop2 <- xtabs( ~. , popdata , subset = pop_frame$v3 == 2 )
nipij_pop1 <- muij_pop1 / sum( muij_pop1 )
nipij_pop2 <- muij_pop2 / sum( muij_pop2 )

# to matrix
attr( muij_pop1 , "call" ) <- NULL
attr( muij_pop2 , "call" ) <- NULL
attr( nipij_pop1 , "call" ) <- NULL
attr( nipij_pop2 , "call" ) <- NULL
muij_pop1 <- unclass( muij_pop1 )
muij_pop2 <- unclass( muij_pop2 )
nipij_pop1 <- unclass( nipij_pop1 )
nipij_pop2 <- unclass( nipij_pop2 )

# build designs
flowdes <-
  sfydesign( ids = ~0 ,
             fpc = ~ 1/prob ,
             data = list( dfa0 , dfa1 ) ,
             nest = TRUE )

# create replicate designs
flowdes_rep <- as.surfrdesign( flowdes , type = "bootstrap" , replicate = 50 )

# estimate flows
gflow_lin_sub1 <- svyflow( ~v0 , subset( flowdes , v1 == 1 ) )
gflow_lin_sub2 <- svyflow( ~v0 , subset( flowdes , v1 == 2 ) )
nflow_lin_sub1 <- svyflow( ~v0 , subset( flowdes , v1 == 1 ) , flow.type = "net" )
nflow_lin_sub2 <- svyflow( ~v0 , subset( flowdes , v1 == 2 ) , flow.type = "net" )
gflow_rep_sub1 <- svyflow( ~v0 , subset( flowdes_rep , v1 == 1 ) )
gflow_rep_sub2 <- svyflow( ~v0 , subset( flowdes_rep , v1 == 2 ) )
nflow_rep_sub1 <- svyflow( ~v0 , subset( flowdes_rep , v1 == 1 ) , flow.type = "net" )
nflow_rep_sub2 <- svyflow( ~v0 , subset( flowdes_rep , v1 == 2 ) , flow.type = "net" )

# test extraction of associated measures
test_that( "extraction of estimates" , {

  # point estimates
  expect_identical( coef( gflow_lin_sub1 ) , surf:::coef.flowstat( gflow_lin_sub1 ) )
  expect_identical( coef( gflow_lin_sub2 ) , surf:::coef.flowstat( gflow_lin_sub2 ) )
  expect_identical( coef( gflow_rep_sub1 ) , surf:::coef.flowstat( gflow_rep_sub1 ) )
  expect_identical( coef( gflow_rep_sub2 ) , surf:::coef.flowstat( gflow_rep_sub2 ) )
  expect_identical( coef( nflow_lin_sub1 ) , surf:::coef.flowstat( nflow_lin_sub1 ) )
  expect_identical( coef( nflow_lin_sub2 ) , surf:::coef.flowstat( nflow_lin_sub2 ) )
  expect_identical( coef( nflow_rep_sub1 ) , surf:::coef.flowstat( nflow_rep_sub1 ) )
  expect_identical( coef( nflow_rep_sub2 ) , surf:::coef.flowstat( nflow_rep_sub2 ) )


  # variances
  expect_identical( vcov( gflow_lin_sub1 ) , surf:::vcov.flowstat( gflow_lin_sub1 ) )
  expect_identical( vcov( gflow_lin_sub2 ) , surf:::vcov.flowstat( gflow_lin_sub2 ) )
  expect_identical( vcov( gflow_rep_sub1 ) , surf:::vcov.flowstat( gflow_rep_sub1 ) )
  expect_identical( vcov( gflow_rep_sub2 ) , surf:::vcov.flowstat( gflow_rep_sub2 ) )
  expect_identical( vcov( nflow_lin_sub1 ) , surf:::vcov.flowstat( nflow_lin_sub1 ) )
  expect_identical( vcov( nflow_lin_sub2 ) , surf:::vcov.flowstat( nflow_lin_sub2 ) )
  expect_identical( vcov( nflow_rep_sub1 ) , surf:::vcov.flowstat( nflow_rep_sub1 ) )
  expect_identical( vcov( nflow_rep_sub2 ) , surf:::vcov.flowstat( nflow_rep_sub2 ) )

  # standard errors
  expect_identical( SE( gflow_lin_sub1 ) , surf:::SE.flowstat( gflow_lin_sub1 ) )
  expect_identical( SE( gflow_lin_sub2 ) , surf:::SE.flowstat( gflow_lin_sub2 ) )
  expect_identical( SE( gflow_rep_sub1 ) , surf:::SE.flowstat( gflow_rep_sub1 ) )
  expect_identical( SE( gflow_rep_sub2 ) , surf:::SE.flowstat( gflow_rep_sub2 ) )
  expect_identical( SE( nflow_lin_sub1 ) , surf:::SE.flowstat( nflow_lin_sub1 ) )
  expect_identical( SE( nflow_lin_sub2 ) , surf:::SE.flowstat( nflow_lin_sub2 ) )
  expect_identical( SE( nflow_rep_sub1 ) , surf:::SE.flowstat( nflow_rep_sub1 ) )
  expect_identical( SE( nflow_rep_sub2 ) , surf:::SE.flowstat( nflow_rep_sub2 ) )

  # coefficients of variation
  expect_identical( cv( gflow_lin_sub1 ) , surf:::cv.flowstat( gflow_lin_sub1 ) )
  expect_identical( cv( gflow_lin_sub2 ) , surf:::cv.flowstat( gflow_lin_sub2 ) )
  expect_identical( cv( gflow_rep_sub1 ) , surf:::cv.flowstat( gflow_rep_sub1 ) )
  expect_identical( cv( gflow_rep_sub2 ) , surf:::cv.flowstat( gflow_rep_sub2 ) )
  expect_identical( cv( nflow_lin_sub1 ) , surf:::cv.flowstat( nflow_lin_sub1 ) )
  expect_identical( cv( nflow_lin_sub2 ) , surf:::cv.flowstat( nflow_lin_sub2 ) )
  expect_identical( cv( nflow_rep_sub1 ) , surf:::cv.flowstat( nflow_rep_sub1 ) )
  expect_identical( cv( nflow_rep_sub2 ) , surf:::cv.flowstat( nflow_rep_sub2 ) )

} )

# test against bias
test_that("compare point estimates vs population values",{
  expect_equivalent( coef( gflow_lin_sub1 ) , muij_pop1 , tolerance = .1 )
  expect_equivalent( coef( gflow_lin_sub2 ) , muij_pop2 , tolerance = .1 )
  expect_equivalent( coef( gflow_rep_sub1 ) , muij_pop1 , tolerance = .1 )
  expect_equivalent( coef( gflow_rep_sub2 ) , muij_pop2 , tolerance = .1 )
  expect_equivalent( coef( nflow_lin_sub1 ) , nipij_pop1 , tolerance = .1 )
  expect_equivalent( coef( nflow_lin_sub2 ) , nipij_pop2 , tolerance = .1 )
  expect_equivalent( coef( nflow_rep_sub1 ) , nipij_pop1 , tolerance = .1 )
  expect_equivalent( coef( nflow_rep_sub2 ) , nipij_pop2 , tolerance = .1 )
} )

# linearized vs replicate
test_that("compare linearized vs replicate: point estimates",{
  expect_identical( coef( gflow_lin_sub1 ) , coef( gflow_rep_sub1 ) )
  expect_identical( coef( gflow_lin_sub2 ) , coef( gflow_rep_sub2 ) )
  expect_identical( coef( nflow_lin_sub1 ) , coef( nflow_rep_sub1 ) )
  expect_identical( coef( nflow_lin_sub2 ) , coef( nflow_rep_sub2 ) )
} )
test_that("compare linearized vs replicate: standard errors",{
  expect_equal( SE( gflow_lin_sub1 ) , SE( gflow_rep_sub1 ) , tolerance = .3 )
  expect_equal( SE( gflow_lin_sub2 ) , SE( gflow_rep_sub2 ) , tolerance = .3 )
  expect_equal( SE( nflow_lin_sub1 ) , SE( nflow_rep_sub1 ) , tolerance = .3 )
  expect_equal( SE( nflow_lin_sub2 ) , SE( nflow_rep_sub2 ) , tolerance = .3 )
} )
test_that("compare linearized vs replicate: coefficients of variation",{
  expect_equal( cv( gflow_lin_sub1 ) , cv( gflow_rep_sub1 ) , tolerance = .3 )
  expect_equal( cv( gflow_lin_sub2 ) , cv( gflow_rep_sub2 ) , tolerance = .3 )
  expect_equal( cv( nflow_lin_sub1 ) , cv( nflow_rep_sub1 ) , tolerance = .3 )
  expect_equal( cv( nflow_lin_sub2 ) , cv( nflow_rep_sub2 ) , tolerance = .3 )
} )
