context("survflow design objects" )

# load data
data( "artificial" )

test_that( "svyflowdesign works on unweighted designs" , {
  flowdes_srs <-
    sfydesign( ids = ~0 ,
               probs = ~ prob ,
               data = list( dfa0 , dfa1 ) ,
               nest = TRUE )
  test_that("summary works on unweighted designs", { summary( flowdes_srs ) } )
  test_that("print works on unweighted designs", { print( flowdes_srs ) } )
} )


test_that( "svyflowdesign works on stratified designs" , {
  flowdes_strat <-
    sfydesign( ids = ~0 ,
               probs = ~ prob ,
               strata = ~ stratum ,
               data = list( surf:::dfa0_strat , surf:::dfa1_strat ) ,
               nest = TRUE )
  test_that("summary works on stratified designs", { summary( flowdes_strat ) } )
  test_that("print works on stratified designs", { print( flowdes_strat ) } )
} )


test_that( "replicate design from unweighted designs" , {
  flowdes_srs <-
    sfydesign( ids = ~0 ,
               probs = ~ prob ,
               data = list( dfa0 , dfa1 ) ,
               nest = TRUE )
  flowdes_srs_rep <- as.surfrdesign( flowdes_srs , type = "bootstrap" , replicate = 50 )
  test_that("summary works on unweighted designs", { summary( flowdes_srs_rep ) } )
  test_that("print works on unweighted designs", { print( flowdes_srs_rep ) } )
} )


test_that( "svyflowdesign works on stratified replicate designs" , {
  flowdes_strat <-
    sfydesign( ids = ~0 ,
               probs = ~ prob ,
               strata = ~ stratum ,
               data = list( surf:::dfa0_strat , surf:::dfa1_strat ) ,
               nest = TRUE )
  flowdes_strat_rep <- as.surfrdesign( flowdes_strat , type = "bootstrap" , replicate = 50 )
  test_that("summary works on stratified replicate designs", { summary( flowdes_strat_rep ) } )
  test_that("summary works on stratified replicate designs", { print( flowdes_strat_rep ) } )
} )

