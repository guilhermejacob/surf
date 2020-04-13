#' @exportS3Method print surflow.design
print.surflow.design <-function(x, varnames=TRUE,design.summaries=FALSE,...){

  if ( "svyrep.design" %in% class(x) ) {
    cat("Call: ")
    print(x$call)
    if (x$type == "Fay")
      cat("Fay's variance method (rho=", x$rho, ") ")
    if (x$type == "BRR")
      cat("Balanced Repeated Replicates ")
    if (x$type == "JK1")
      cat("Unstratified cluster jacknife (JK1) ")
    if (x$type == "JKn")
      cat("Stratified cluster jackknife (JKn) ")
    if (x$type == "bootstrap")
      cat("Survey bootstrap ")
    if (x$type == "mrbbootstrap")
      cat("Multistage rescaled bootstrap ")
    if (x$type == "subbootstrap")
      cat("(n-1) bootstrap ")
    nweights <- ncol(x$repweights)
    cat("with", nweights, "replicates")
    if (!is.null(x$mse) && x$mse)
      cat(" and MSE variances")
    cat(".\n")
  } else {
    n<-NROW(x$cluster)
    if (x$has.strata) cat("Stratified ")
    un<-length(unique(x$cluster[,1]))
    if(n==un){
      cat("Independent Sampling design")
      is.independent<-TRUE
      if (is.null(x$fpc$popsize))
        cat(" (with replacement)\n")
      else cat("\n")
    } else {
      cat(NCOL(x$cluster),"- level Cluster Sampling design")
      if (is.null(x$fpc$popsize))
        cat(" (with replacement)\n")
      else cat("\n")
      nn<-lapply(x$cluster,function(i) length(unique(i)))
      cat(paste("With (",paste(unlist(nn),collapse=", "),") clusters.\n",sep=""))
      is.independent<-FALSE
    }

    print(x$call)
    if (design.summaries){
      cat("Probabilities:\n")
      print(summary(x$prob))
      if(x$has.strata){
        if (NCOL(x$cluster)>1)
          cat("First-level ")
        cat("Stratum Sizes: \n")
        oo<-order(unique(x$strata[,1]))
        a<-rbind(obs=table(x$strata[,1]),
                 design.PSU=x$fpc$sampsize[!duplicated(x$strata[,1]),1][oo],
                 actual.PSU=table(x$strata[!duplicated(x$cluster[,1]),1]))
        print(a)
      }
      if (!is.null(x$fpc$popsize)){
        if (x$has.strata) {
          cat("Population stratum sizes (PSUs): \n")
          s<-!duplicated(x$strata[,1])
          a<-x$fpc$popsize[s,1]
          names(a)<-x$strata[s,1]
          a<-a[order(names(a))]
          print(a)
        } else {
          cat("Population size (PSUs):",x$fpc$popsize[1,1],"\n")
        }
      }
    }
  }
  cat("Number of repetitions:", length( x$variables ) - 1 , "\n" )
  if (varnames){
    cat("Data variables:\n")
    print( Reduce( intersect , sapply( x$variables , colnames ) ) )
  }
  invisible(x)
}

#' @exportS3Method summary surflow.design
summary.surflow.design <-function(object,varnames=TRUE,...){

  if ( "svyrep.design" %in% class(object) ) {
    cat("Call: ")
    print(object$call)
    if (object$type == "Fay")
      cat("Fay's variance method (rho=", object$rho, ") ")
    if (object$type == "BRR")
      cat("Balanced Repeated Replicates ")
    if (object$type == "JK1")
      cat("Unstratified cluster jacknife (JK1) ")
    if (object$type == "JKn")
      cat("Stratified cluster jackknife (JKn) ")
    if (object$type == "bootstrap")
      cat("Survey bootstrap ")
    if (object$type == "mrbbootstrap")
      cat("Multistage rescaled bootstrap ")
    if (object$type == "subbootstrap")
      cat("(n-1) bootstrap ")
    nweights <- ncol(object$repweights)
    cat("with", nweights, "replicates")
    if (!is.null(object$mse) && object$mse)
      cat(" and MSE variances")
    cat(".\n")
  } else {
    n<-NROW(object$cluster)
    if (object$has.strata) cat("Stratified ")
    un<-length(unique(object$cluster[,1]))
    if(n==un){
      cat("Independent Sampling design")
      is.independent<-TRUE
      if (is.null(object$fpc$popsize))
        cat(" (with replacement)\n")
      else cat("\n")
    } else {
      cat(NCOL(object$cluster),"- level Cluster Sampling design")
      if (is.null(object$fpc$popsize))
        cat(" (with replacement)\n")
      else cat("\n")
      nn<-lapply(object$cluster,function(i) length(unique(i)))
      cat(paste("With (",paste(unlist(nn),collapse=", "),") clusters.\n",sep=""))
      is.independent<-FALSE
    }

    print(object$call)
    if (TRUE){
      cat("Probabilities:\n")
      print(summary(object$prob))
      if(object$has.strata){
        if (NCOL(object$cluster)>1)
          cat("First-level ")
        cat("Stratum Sizes: \n")
        oo<-order(unique(object$strata[,1]))
        a<-rbind(obs=table(object$strata[,1]),
                 design.PSU=object$fpc$sampsize[!duplicated(object$strata[,1]),1][oo],
                 actual.PSU=table(object$strata[!duplicated(object$cluster[,1]),1]))
        print(a)
      }
      if (!is.null(object$fpc$popsize)){
        if (object$has.strata) {
          cat("Population stratum sizes (PSUs): \n")
          s<-!duplicated(object$strata[,1])
          a<-object$fpc$popsize[s,1]
          names(a)<-object$strata[s,1]
          a<-a[order(names(a))]
          print(a)
        } else {
          cat("Population size (PSUs):",object$fpc$popsize[1,1],"\n")
        }
      }
    }
  }
  cat("Number of repetitions:", length( object$variables ) - 1 , "\n" )
  if (varnames){
    cat("Data variables:\n")
    print( Reduce( intersect , sapply( object$variables , colnames ) ) )
  }
  invisible(object)
}

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

  # format output
  cmat <- format( cmat[,] ,
                  digits = ifelse( attr( x , "statistic" ) == "muij" , 0 , 4 ) ,
                  nsmall = ifelse( attr( x , "statistic" ) == "muij" , 0 , 4 ) ,
                  scientific = FALSE ,
                  trim = TRUE ,
                  justify = "centre" , drop0trailing = TRUE )
  semat <- format( semat[,] ,
                   digits = ifelse( attr( x , "statistic" ) == "muij" , 2 , 4 ) ,
                   nsmall = ifelse( attr( x , "statistic" ) == "muij" , 2 , 4 ) ,
                   scientific = FALSE ,
                   trim = TRUE ,
                   justify = "centre" , drop0trailing = TRUE )
  cvmat <- format( cvmat[,] ,
                   digits = 6 ,
                   nsmall = 6 ,
                   scientific = FALSE ,
                   trim = TRUE ,
                   justify = "centre" , drop0trailing = TRUE )

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
  print( oomat , quote = FALSE )

}

#' @exportS3Method print flowstat
print.flowstat <- function( x , var.type = c("se","var","cv") , ... ) {

  cat( paste0( "Model " , x$model , "\n" ) )
  cat( paste0( "\nInitial Response Probability" , "\n" ) )
  print( x[["psi"]] )
  cat( paste0( "\nRespondent to Respondent Transition Probability" , "\n" ) )
  print( x[["rhoRR"]] )
  cat( paste0( "\nNon-Respondent to Non-Respondent Transition Probability" , "\n" ) )
  print( x[["rhoMM"]] )
  cat( paste0( "\nGross Flows" , "\n" ) )
  print.svymstat( x[["muij"]] )
  invisible(x)

}
