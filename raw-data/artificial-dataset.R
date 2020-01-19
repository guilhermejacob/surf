# código for simple random sample without replacement

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

# adds non-response
pop_frame[ as.logical( rbinom(N,1,.2) ) , 1 ] <- NA
pop_frame[ !is.na( pop_frame[,1] ) & as.logical( rbinom(N,1,.1) ) , 2 ] <- NA
pop_frame[ is.na( pop_frame[,2] ) & as.logical( rbinom(N,1,.7) ) , 1 ] <- NA
pop_frame <- apply(pop_frame,2,factor)

# extract sample
smp_frame <- pop_frame[ as.logical( sampling::srswor( n , nrow( pop_frame ) ) ) , ]
rm( pop_frame )

# separa em dois
dfa0 <- data.frame( v0 = smp_frame[ , 1 ] )
dfa0$prob <- n/N
dfa1 <- data.frame( v0 = smp_frame[ , 2 ] )

# salva resultado
save( list = c("dfa0","dfa1") , file = "data/artificial.rda" )
