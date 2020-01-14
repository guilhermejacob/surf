## print methods

print.survflow.design<-function(x,varnames=FALSE,design.summaries=FALSE,...){
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
  cat("Number of repetitions:", length( attr( x , "data.pairs" ) ) , "\n" )
  if (varnames){
    cat("Data variables:\n")
    Reduce( intersect , c( x$variables , sapply( attr( x , "data.pairs" ), colnames ) ) )
  }
  invisible(x)
}
