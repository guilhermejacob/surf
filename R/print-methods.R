#' @export
#' @method print survflow.design
print.survflow.design <-function(x,varnames=TRUE,design.summaries=FALSE,...){

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

#' @export
#' @method print flowstat
print.flowstat <- function( object , digits = 0 , ...) {
  if ( attr( object , "statistic" ) == "gross" ) { object <- round( object , digits = digits ) ; attr( object , "var" ) <- round( attr( object , "var" ) , digits = digits ) }
  # print( object[,] )
  printCoefmat( object )
}
