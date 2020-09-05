# create matrix of counts
LFS79.0809.counts <- matrix( c(9222L, 221L, 256L, 996L, 128L, 322L, 164L, 69L, 662L, 151L, 5941L, 676L, 473L, 59L, 292L, 4353L) , ncol = 4 , byrow = FALSE )
LFS79.0910.counts <- matrix( c(9697L, 177L, 326L, 554L, 169L, 317L, 159L, 59L, 355L, 143L, 6522L, 339L, 474L, 46L, 423L, 4225L) , ncol = 4 , byrow = FALSE )

# set classifications
colnames( LFS79.0809.counts ) <- rownames( LFS79.0809.counts ) <- c( "E" , "U" , "N" , NA )
colnames( LFS79.0910.counts ) <- rownames( LFS79.0910.counts ) <- c( "E" , "U" , "N" , NA )

# transform into synthetic sample
df.list <-
  lapply( list( LFS79.0809.counts , LFS79.0910.counts ) , function ( cmat ) {

  # set up data.frame
  ll <- expand.grid( i = c( "E" , "U" , "N" , NA ) , j = c( "E" , "U" , "N" , NA ) )
  ll$ijcode <- seq_len( nrow(ll) )
  ll$counts <- as.vector( cmat )

  df.art <- data.frame( ijcode = unlist( lapply( seq_len( nrow( ll ) ) , function( k ) rep( ll$ijcode[k] , ll$counts[k] ) ) ) )
  df.art <- merge( df.art , ll[ , 1:3] , by = "ijcode" , sort = FALSE , all.x = TRUE )

  # set random order
  df.art <- df.art[ sample.int( nrow( df.art ) , nrow(df.art ) , replace = FALSE ) , ]

  # remove rownames
  rownames( df.art ) <- NULL

  # reformat data.frame
  df.art$ijcode <- NULL
  colnames( df.art ) <- c( "y1" , "y2" )

  # add weights
  df.art$prob <- 1

  # return final data.frame
  df.art

} )

# test counts
table( df.list[[1]]$y1 , df.list[[1]]$y2 , useNA = "always" )
table( df.list[[2]]$y1 , df.list[[2]]$y2 , useNA = "always" )

# save data
LFS79.0809 <- df.list[[1]]
save( list = "LFS79.0809" , file = "data/LFS79.0809.rda" )
LFS79.0910 <- df.list[[2]]
save( list = "LFS79.0910" , file = "data/LFS79.0910.rda" )
