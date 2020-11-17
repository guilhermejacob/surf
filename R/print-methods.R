#' @exportS3Method print svymstat
print.svymstat <- function( x , var.type = c("se","var","cv") , ... ) {

  # variance type
  var.type <- match.arg( var.type , several.ok = FALSE )

  # collect coefficients
  cmat <- x
  attr( cmat , "var" ) <- NULL
  attr( cmat , "statistic" ) <- NULL
  attr( cmat , "class" ) <- NULL


  # collect standard-errors, variance, and cv
  vmat <- attr( x , "var" )
  semat <- sqrt( vmat )
  cvmat <- semat / cmat

  # rounding
  vmat <- round( vmat , ifelse( attr( x , "statistic" ) == "muij" , 0 , 4 ) )
  semat <- round( semat , ifelse( attr( x , "statistic" ) == "muij" , 0 , 4 ) )
  cvmat <- round( cvmat , ifelse( attr( x , "statistic" ) == "muij" , 0 , 4 ) )
  cmat <- round( cmat , ifelse( attr( x , "statistic" ) == "muij" , 0 , 4 ) )

  # format output
  cmat <- format( cmat[,] ,
                  digits = ifelse( attr( x , "statistic" ) == "muij" , 2 , 4 ) ,
                  nsmall = ifelse( attr( x , "statistic" ) == "muij" , 2 , 4 ) ,
                  scientific = FALSE ,
                  trim = FALSE ,
                  justify = "centre" , drop0trailing = FALSE )
  semat <- format( semat[,] ,
                   digits = ifelse( attr( x , "statistic" ) == "muij" , 2 , 4 ) ,
                   nsmall = ifelse( attr( x , "statistic" ) == "muij" , 2 , 4 ) ,
                   scientific = FALSE ,
                   trim = FALSE ,
                   justify = "centre" , drop0trailing = TRUE )
  cvmat <- format( cvmat[,] ,
                   digits = 2 ,
                   nsmall = 2 ,
                   scientific = FALSE ,
                   trim = FALSE ,
                   justify = "centre" , drop0trailing = FALSE )

  # get header
  if ( attr( x , "statistic" ) == "muij" ) opheader <- "gross flows" else if ( attr( x , "statistic" ) == "pij" ) opheader <- "transition probabilities"

  # # print flow estimates header
  # cat(paste0(opheader , "\nestimates\n" ) )

  # print estimates
  print( cmat , quote = FALSE )

  # print standard-errors
  cat(paste0("\n",
             switch( var.type ,
                     se = "SE" ,
                     var = "variances" ,
                     cv = "coefficients of variation" )
             , "\n"))

  oomat <- switch( var.type ,
                   se = semat ,
                   var = var ,
                   cv = cvmat )
  print( oomat , quote = FALSE , digits = 4 )

}

#' @exportS3Method print flowstat
print.flowstat <- function( x , ... ) {

  # print model type
  cat( paste0( "Model " , x$model , "\n" ) )

  # print estimates of non-response mechanism
  if ( !is.null( x[["psi"]] ) ){
    cat( paste0( "\nInitial Response Probability" , "\n" ) )
    print( x[["psi"]] )
    cat( paste0( "\nRespondent to Respondent Transition Probability" , "\n" ) )
    print( x[["rho"]] )
    cat( paste0( "\nNon-Respondent to Non-Respondent Transition Probability" , "\n" ) )
    print( x[["tau"]] )
  }

  # print eta
  # cat( paste0( "\nInitial Distribution" , "\n" ) )
  # print( x[["eta"]] )

  # print gamma
  # cat( paste0( "\nFinal Distribution" , "\n" ) )
  # print( x[["gamma"]] )

  # print gross flows
  cat( paste0( "\nGross Flows" , "\n" ) )
  print.svymstat( x[["muij"]] )

  # print model fit
  if ( !is.null( attr( x , "adj.chisq" ) ) ){
    cat( paste0( "\n" , "\n" ) )
    print( attr( x , "adj.chisq" ) )
  }

  # return
  invisible(x)

}
