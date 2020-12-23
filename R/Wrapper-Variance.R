variance_fun <- function( xx , ww , res , design ) {

  switch( res$model ,
          A = { variance.fun <- modelA.variance } ,
          B = { variance.fun <- modelB.variance } ,
          C = { variance.fun <- modelC.variance } ,
          D = { variance.fun <- modelD.variance } )

  variance.fun( xx , ww , res , design )

}

