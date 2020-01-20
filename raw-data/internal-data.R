# code for stratified sample

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

# create strata
pop_frame$strat <- ifelse( pop_frame$v0 == 1 , rbinom( nrow( pop_frame ) , 1 , .7 ) , rbinom( nrow( pop_frame ) , 1 , .3 ) ) + 1
pop_frame <- pop_frame[ order( pop_frame$strat ) , ]
pop_frame <- data.frame( apply(pop_frame[,] , 2 , factor) , stringsAsFactors = TRUE )

# adds non-response
pop_frame[ as.logical( rbinom(N,1,.2) ) , 1 ] <- NA
pop_frame[ !is.na( pop_frame[,1] ) & as.logical( rbinom(N,1,.1) ) , 2 ] <- NA
pop_frame[ is.na( pop_frame[,2] ) & as.logical( rbinom(N,1,.7) ) , 1 ] <- NA

# extract sample
smp_info <- sampling::strata( pop_frame , stratanames = "strat" , size = round( c( .05 , .10 ) * table( pop_frame[,"strat"] ) ) , method = "srswor" )
smp_frame <- cbind( pop_frame[ smp_info$ID_unit , 1:2 ] , smp_info[ , -1:-2 ] )
rm( pop_frame )

# adjusts column names
colnames( smp_frame ) <- tolower( colnames( smp_frame ) )

# splits data.frame
dfa0_strat <- smp_frame[ , -2 ]
dfa1_strat <- smp_frame[ , 2 , drop = FALSE ]
colnames( dfa1_strat ) <- "v0"

# save results
devtools::use_data( muij_pop , nipij_pop , dfa0_strat , dfa1_strat , internal = TRUE , overwrite = TRUE , compress = "xz" )
