# código para pseudo população

# gerador aleatório
set.seed(123)

### PARÂMETROS DA SIMULAÇÃO

# define parâmetros da população artificial
N = 10^5
eta_pop <- c( .9 , .05, .05 )
pij_pop <- matrix( c(.80, .15, .05, .30, .60, .10, .10, .10, .80 ) , ncol = 3 , nrow = 3 , byrow = T )
nipij_pop <- sweep( pij_pop , 1 , eta_pop , "*" )
muij_pop <- nipij_pop * N

# matriz de classificações
class_table <- expand.grid( data.frame( v0 = seq_len( nrow(nipij_pop) ) , v1 = seq_len( nrow(nipij_pop) ) ) )
class_table <- class_table[ order( class_table$v0 ) , ]
class_table[ ,"k_ij" ] <- as.character( seq_len( nrow( class_table ) ) )

# define tamaho da amostra
n = as.integer(10^4)

# padrão de não resposta
pop_psi <- .8
pop_rhoRR <- .9
pop_rhoMM <- .7

### GERA POPULAÇÃO ARTIFICIAL

# gera população de tamanho N
pop_fullresponse <- t( rmultinom( as.integer( N ) , size = 1 , prob = as.numeric( t( nipij_pop ) ) ) )

# transforma em tipologia de transições
pop_fullresponse <- apply( pop_fullresponse , 1 , function( z ) seq_len( ncol( pop_fullresponse ) )[ as.logical( z ) ] )

# transforma em transição de categorias
pop_fullresponse <- data.frame( "id" = seq_len( N ) , "k_ij" = pop_fullresponse , row.names = NULL , stringsAsFactors = FALSE )

# combina com classificação de transições
pop_fullresponse <- merge( pop_fullresponse , class_table , by = c( "k_ij" ) , all.x = TRUE , all.y = FALSE , sort = TRUE )

# reordena para posoções originais
pop_fullresponse <- pop_fullresponse[ order( pop_fullresponse$id ) , ]

# remove nome das linhas
rownames( pop_fullresponse ) <- NULL

# transforma categorias em factors
pop_fullresponse[, c( "v0" , "v1" ) ] <- lapply( pop_fullresponse[, c( "v0" , "v1" ) ] , factor , levels = c( 1:3 ) , labels = 1:3 )

### APLICA PROCESSO ESTOCÁSTICO DE NÃO-RESPOSTA

# copia população com resposta completa
pop_nonrespose <- pop_fullresponse

# adiciona não-resposta aleatória no tempo t-1
pop_nonrespose[ as.logical( rbinom( nrow( pop_fullresponse ) , 1 , 1 - pop_psi ) ) , "v0" ] <- NA

# adiciona não-resposta no tempo t para os indivíduos respondentes no tempo t-1
pop_nonrespose[ !is.na( pop_nonrespose[,"v0"] ) & as.logical( rbinom( nrow( pop_nonrespose ) , 1 , 1 - pop_rhoRR ) ) , "v1" ] <- NA

# adiciona não-resposta no tempo t para os indivíduos não-respondentes no tempo t-1
pop_nonrespose[  is.na( pop_nonrespose[,"v0"] ) & as.logical( rbinom( nrow( pop_nonrespose ) , 1, pop_rhoMM ) ) , "v1" ] <- NA

### EXTRAI AMOSTRA

# extrai amostra de n observações
smp_df <- pop_nonrespose[ as.logical( rbinom( n , 1 , n / N ) ) , c( "v0" , "v1" ) ]

# cria informações do desenho amostral
smp_df$prob <- N / n            # peso
smp_df$fpcs <- as.integer( N )  # correção de população finita

# CRIA INFORMAÇÕES

# cria bases de dados pareadas
dfa0 <- smp_df[ , -2 , drop = FALSE ]
dfa1 <- smp_df[ , 2 , drop = FALSE ]
colnames( dfa0 )[1] <- "y"
colnames( dfa1 )[1] <- "y"

# remove objetos
rm( pop_fullresponse , pop_nonrespose , smp_df ) ; gc()

# save results
save( list = c("dfa0","dfa1") , file = "data/artificial.rda" , compress = "xz" )
