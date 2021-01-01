linearization_fun <- function( xx , ww , res , design ) {

  switch( res$model ,
          A = { linearization.fun <- modelA.linearization } ,
          B = { linearization.fun <- modelB.linearization } ,
          C = { linearization.fun <- modelC.linearization } ,
          D = { linearization.fun <- modelD.linearization } )

  linearization.fun( xx , ww , res , design )

}
