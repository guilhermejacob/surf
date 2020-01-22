context("survflow design objects" )

# define seed for random numbers
set.seed(123)

# load data
data( "artificial" )

# create designs
flowdes_srs <-
  sfydesign( ids = ~0 ,
             data = list( dfa0 , dfa1 ) ,
             nest = TRUE )
flowdes_strat <-
  sfydesign( ids = ~0 ,
             probs = ~ prob ,
             strata = ~ stratum ,
             data = list( surf:::dfa0_strat , surf:::dfa1_strat ) ,
             nest = TRUE )
flowdes_srs_rep <- as.surfrdesign( flowdes_srs , type = "bootstrap" , replicate = 50 )
flowdes_strat_rep <- as.surfrdesign( flowdes_strat , type = "bootstrap" , replicate = 50 )

test_that( "svyflowdesign works on unweighted designs" , {
  verify_output("summary_unweighted_lin.txt", { summary( flowdes_srs ) } )
  verify_output("print_unweighted_lin.txt", { print( flowdes_srs ) } )
} )


test_that( "svyflowdesign works on stratified designs" , {
  verify_output("summary_weighted_lin.txt", { summary( flowdes_strat ) } )
  verify_output("print_weighted_lin.txt", { print( flowdes_strat ) } )
} )


test_that( "replicate design from unweighted designs" , {
  verify_output("summary_unweighted_rep.txt", { summary( flowdes_srs_rep ) } )
  verify_output("print_unweighted_rep.txt", { print( flowdes_srs_rep ) } )
} )


test_that( "svyflowdesign works on stratified replicate designs" , {
  verify_output("summary_weighted_rep.txt", { summary( flowdes_strat_rep ) } )
  verify_output("print_weighted_rep.txt", { print( flowdes_strat_rep ) } )
} )

