#' @name svyflow
#' @title Gross flow estimation between categories
#'
#' @description Compute gross flows for data from complex surveys with repeated samples.
#'
#' @param x  a one-sided formula indicating a \emph{factor} variable.
#' @param design  survey design object
#' @param model  model for non-response. Possibilities: \code{"A", "B", "C", "D"}. Defaults to \code{model = "A"}.
#' @param tol  Tolerance for iterative proportional fitting. Defaults to \code{1e-4}.
#' @param maxit  Maximum number of iterations for iterative proportional fitting. Defaults to \code{maxit = 5000}.
#' @param verbose  Print proportional fitting iterations. Defaults to \code{verbose = FALSE}.
#' @param as.zero.flows  Should zeroes in the observed gross flows should be considered as zeroes in the population transition probability matrix? Defaults to \code{as.zero.flows = FALSE}.
#' @param influence  Should influence functions estimates be stored? Defaults to \code{influence = FALSE}.
#' @param ...  future expansion.
#'
#' @details It is important to distinguish "missing" responses from "unnaplicable" responses. This is feasible by subsetting the design
#' for only applicable responses (with actual missing responses, if that is the case). For instance, suppose that we have two variables encoded as
#' employed/unemployed, with NAs if the response is missing or is unnaplicable. An \code{NA} might be a person that did not respond \emph{or} a person
#' who was under the working-age at the time of the survey. It is important to distinguish across those, as only one of those cases is an
#' actual non-response. You could do that by looking for people who were in working age in any round, for instance. This can be done by using \code{subset},
#' as you should for a \code{survey design} object.
#'
#' @return Objects of class "flowstat", which are tables with a "var" attribute giving the variance and a "statistic" attribute giving the type of flow.
#'
#' These objects have methods for coef, vcov, SE, and cv.
#'
#' @author Guilherme Jacob
#'
#' @examples
#'
#' @references GUTIERREZ, H. A.; TRUJILLO, L.; SILVA, P. L. N. The estimation of gross flows in complex surveys with random nonresponse.
#' \emph{Survey Methodology}, v. 40, n. 2, p. 285–321, dec. 2014. URL \url{https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X201400214113}.
#'
#' LUMLEY, T. \emph{Complex Surveys:} A guide to analysis using R.
#' Hoboken: John Wiley & Sons, 2010. (Wiley Series in Survey Methodology). ISBN 978-0-470-28430-8.
#'
#' @examples
#'
#' # load library
#' library( survey )
#' library( surf )
#'
#' # load data
#' data( "LFS79.0809" )
#'
#' # create surf design object
#' lfs.des <- svydesign( ids = ~0 , probs = ~ prob , data = LFS79.0809 , nest = TRUE )
#'
#' # flow estimates
#' estflows <- svyflow( ~y1+y2 , design = lfs.des )
#' coef( estflows$muij )
#' SE( estflows$muij )
#'
#' @export
#' @rdname svyflow
#' @method svyflow survey.design2
svyflow.survey.design2 <- function( x , design , model = c("A","B","C","D") , tol = 1e-4 , maxit = 5000 , verbose = FALSE , as.zero.flows = FALSE , influence = FALSE , ... ){

  # test values
  model <- match.arg( model , several.ok = FALSE )

  # collect sample data and put in single data frame
  xx <- stats::model.frame( x , data = design$variables , na.action = stats::na.pass )

  # test column format
  if ( !all( sapply( xx , is.factor ) ) ) stop( "this function is only valid for factors." )

  # Gets levels of factors for both time periods
  xlevels <- lapply( xx , function(zz) levels( zz ) )

  # gets dimension of variable for which flows are to be estimated
  ll <- lapply( xlevels , function(zz) length( zz ) )
  if ( !all( ll[[1]] == ll[[2]] & xlevels[[1]] == xlevels[[2]] ) ) stop( "inconsistent categories across rounds." )

  # when levels are the same across periods, returns unique levels of variable for which flows are to be estimated
  xlevels <- unique( unlist(xlevels) )

  # check for ordered categories
  ll <- unique( unlist(ll) )
  has.order <- ifelse( all( sapply( xx , is.ordered ) ) , TRUE , FALSE )

  # collect weights
  ww <-  stats::weights( design )

  # estimate counts
  Amat <- stats::xtabs( c(ww,0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE , drop.unused.levels = FALSE )

  # non-response cells
  Nij <- Amat[ -nrow( Amat ) , -ncol( Amat ) ]
  Ri <- Amat[ , ncol(Amat) ][ - nrow( Amat ) ]
  Cj <- Amat[ nrow( Amat ) , ][ - ncol( Amat ) ]
  M <- Amat[ nrow( Amat ) , ncol( Amat ) ]

  # test for zero counts in observed flows
  if ( any( Nij <= 0 ) & !as.zero.flows ) {
    print( Amat , quote = FALSE )
    stop( "Some observed flow cells had zero counts.
             If those are zero counts in the population transition probability matrix, consider using as.zero.flows = TRUE.
             If not, consider collapsing categories." )
  } else if ( any( Nij <= 0 ) & as.zero.flows ) {
    warning( "Some observed flow cells had zero counts. Model fitted with zeroes in the population transition probability matrix." , immediate. = TRUE )
  }

  # treat full response
  if ( all( c( Ri , Cj , M ) == 0 ) ) {

    # model fitting
    N <- sum( Nij )
    mfit <-
      list( "observed.counts" = Amat ,
            eta = rowSums( Nij ) / N ,
            pij = sweep( Nij , 1 , rowSums( Nij ) , "/" ) )
    mfit$nipij <- sweep( mfit$pij , 1 , mfit$eta , "*" )
    mfit$gamma <- colSums( mfit$nipij )
    mfit$delta <- N * ( mfit$gamma - mfit$eta )
    mfit$muij <- N*sweep( mfit$pij , 1 , mfit$eta , "*" )

    # visit proportion
    model.expected <- rbind( cbind( mfit$nipij , rep( 0 , ncol( Nij ) ) ) , rep( 0 , nrow( Nij ) + 1 ) )
    interview.number <- data.frame( visit = rep( 1 , length( ww ) ) )
    vVec <- c( stats::xtabs( ww ~ visit , data = interview.number , addNA = TRUE , drop.unused.levels = FALSE ) / sum( ww ) )
    model.expected <- outer( model.expected , vVec )

    # rao-scott adjustment of the chi-square
    K <- dim( model.expected )[1] - 1
    n.parms <- K + K^2
    pearson <- rao.scott( xx , interview.number , ww , model.expected , design , n.parms , pij.zero = sum( mfit$pij == 0 ) )

    # linearization
    llin <- FullResponse.linearization( xx , ww , mfit , design )

    # calculate variance
    mvar <- lapply( llin , function(z) survey::svyrecvar( sweep( z , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )

    # build results list
    res <- sapply( c( "eta" , "pij" , "muij" , "gamma" , "delta" ) , function(z) {
      if ( z %in% c( "psi" , "rho" , "tau" , "eta" , "gamma" , "delta" ) ) {
        this.stats <- mfit[[z]]
        attr( this.stats , "var" ) <- mvar[[z]]
        names( this.stats ) <- if ( length( this.stats ) > 1 ) xlevels else z
        colnames( attr( this.stats , "var" ) ) <- rownames( attr( this.stats , "var" ) ) <- if ( length( attr( this.stats , "var" ) ) > 1 ) xlevels else z
        class( this.stats ) <- "svystat"
        attr( this.stats , "statistic" ) <- z
      } else if ( z %in% c( "pij" , "muij" ) ) {
        this.stats <- mfit[[z]]
        these.classes <- expand.grid( dimnames( this.stats ) )
        these.classes <- these.classes[ order( these.classes[, 2 ] ) , ]
        these.classes <- apply( these.classes , 1 , paste, collapse = ":" )
        this.vector <- c( t( this.stats ) )
        names( this.vector ) <- these.classes
        this.vmat <-  mvar[[z]]
        dimnames( this.vmat ) <- list( these.classes , these.classes )
        attr( this.vector , "var" ) <- this.vmat
        attr( this.vector , "categories" ) <- dimnames( this.stats )
        class( this.vector ) <- c( "svymstat" , "svystat" )
        attr( this.vector , "statistic" ) <- z
        this.stats <- this.vector
      }
      return(this.stats)
    } , simplify = FALSE )

    # add influence attribute
    for ( this.stat in names( res ) ) attr( res[[ this.stat ]] , "influence" ) <- if ( influence ) llin[[this.stat]] else NULL

    # create final object
    rval <- res[ c( "eta" , "gamma" , "pij" , "muij" , "delta" ) ]
    rval$model <- "Full Response"
    class(rval) <- "flowstat"
    attr( rval , "formula" ) <- x
    attr( rval , "has.order" ) <- has.order
    attr( rval , "iter" ) <- NA
    attr( rval , "adj.chisq" ) <- pearson

    # return final object
    return( rval )

  }

  # model fitting
  mfit <- ipf( Amat , model , tol = tol , maxit = maxit , verbose = verbose , keep.info = FALSE )

  # visit proportions
  interview.number <- data.frame( visit = rep( 1 , length( ww) ) )
  vVec <- c( stats::xtabs( ww ~ visit , data = interview.number , addNA = TRUE , drop.unused.levels = FALSE ) / sum( ww ) )
  model.expected <- outer( mfit$estimated.props , vVec )

  # rao-scott adjustment of the chi-square
  K <- dim( model.expected )[1] - 1
  n.parms <- switch (model, A = { K^2 + K + 3 } , B = { K^2 + 2*K + 2 } , C = { K^2 + 3*K + 1 } , D = { K^2 + 3*K + 1 } )
  pearson <- rao.scott( xx , interview.number , ww , model.expected , design , n.parms , pij.zero = sum( mfit$pij == 0 ) )

  # estimate linearized variables
  llin <- linearization_fun( xx , ww , res = mfit , design = design )
  llin[ lapply( llin , class ) %in% "numeric" ] <- lapply( llin[ lapply( llin , class ) %in% "numeric" ] , matrix , nrow = length( ww ) )

  # calculate variance
  mvar <- lapply( llin , function(z) survey::svyrecvar( sweep( z , 1 , ww , "*" ) , clusters = design$cluster , stratas = design$strata , fpcs = design$fpc , postStrata = design$postStrata ) )

  # build results list
  res <- sapply( c( "psi" , "rho" , "tau" , "eta" , "pij" , "muij" , "gamma" , "delta" ) , function(z) {
    if ( z %in% c( "psi" , "rho" , "tau" , "eta" , "gamma" , "delta" ) ) {
      this.stats <- mfit[[z]]
      attr( this.stats , "var" ) <- mvar[[z]]
      names( this.stats ) <- if ( length( this.stats ) > 1 ) xlevels else z
      colnames( attr( this.stats , "var" ) ) <- rownames( attr( this.stats , "var" ) ) <- if ( length( attr( this.stats , "var" ) ) > 1 ) xlevels else z
      class( this.stats ) <- "svystat"
      attr( this.stats , "statistic" ) <- z
    } else if ( z %in% c( "pij" , "muij" ) ) {
      this.stats <- mfit[[z]]
      these.classes <- expand.grid( dimnames( this.stats ) )
      these.classes <- these.classes[ order( these.classes[, 2 ] ) , ]
      these.classes <- apply( these.classes , 1 , paste, collapse = ":" )
      this.vector <- c( t( this.stats ) )
      names( this.vector ) <- these.classes
      this.vmat <-  mvar[[z]]
      dimnames( this.vmat ) <- list( these.classes , these.classes )
      attr( this.vector , "var" ) <- this.vmat
      attr( this.vector , "categories" ) <- dimnames( this.stats )
      class( this.vector ) <- c( "svymstat" , "svystat" )
      attr( this.vector , "statistic" ) <- z
      this.stats <- this.vector
    }
    return(this.stats)
  } , simplify = FALSE )

  # add influence attribute
  for ( this.stat in names( res ) ) attr( res[[ this.stat ]] , "influence" ) <- if ( influence ) llin[[this.stat]] else NULL

  # create final object
  rval <- res[ c( "psi" , "rho" , "tau" , "eta" , "gamma" , "pij" , "muij" , "delta" ) ]
  rval$model <- mfit$model
  class(rval) <- "flowstat"
  attr( rval , "formula" ) <- x
  attr( rval , "has.order" ) <- has.order
  attr( rval , "iter" ) <- mfit$iter
  attr( rval , "adj.chisq" ) <- pearson

  # return final object
  rval

}

#' @export
#' @rdname svyflow
#' @method svyflow svyrep.design
svyflow.svyrep.design <- function( x , design , model = c("A","B","C","D") , tol = 1e-4 , maxit = 5000 , verbose = FALSE , as.zero.flows = FALSE , influence = FALSE , ... ){

  # return warning
  warning( "svyflow is experimental for replicate-based survey design objects." , immediate. = TRUE )

  # test values
  model <- match.arg( model , several.ok = FALSE )

  # collect sample data and put in single data frame
  xx <- stats::model.frame( x , data = design$variables , na.action = stats::na.pass )

  # test column format
  if ( !all( sapply( xx , is.factor ) ) ) stop( "this function is only valid for factors." )

  # Gets levels of factors for both time periods
  xlevels <- lapply( xx , function(zz) levels( zz ) )

  # gets dimension of variable for which flows are to be estimated
  ll <- lapply( xlevels , function(zz) length( zz ) )
  if ( !all( ll[[1]] == ll[[2]] & xlevels[[1]] == xlevels[[2]] ) ) stop( "inconsistent categories across rounds." )

  # when levels are the same across periods, returns unique levels of variable for which flows are to be estimated
  xlevels <- unique( unlist(xlevels) )

  # check for ordered categories
  ll <- unique( unlist(ll) )
  has.order <- ifelse( all( sapply( xx , is.ordered ) ) , TRUE , FALSE )

  # collect weights
  ww <-  stats::weights( design , "sampling" )

  # estimate counts
  Amat <- stats::xtabs( c(ww,0) ~ . , data = rbind(xx , rep(NA,ncol(xx))) , addNA = TRUE , drop.unused.levels = FALSE )

  # non-response cells
  Nij <- Amat[ -nrow( Amat ) , -ncol( Amat ) ]
  Ri <- Amat[ , ncol(Amat) ][ - nrow( Amat ) ]
  Cj <- Amat[ nrow( Amat ) , ][ - ncol( Amat ) ]
  M <- Amat[ nrow( Amat ) , ncol( Amat ) ]

  # test for zero counts in observed flows
  if ( any( Nij <= 0 ) & !as.zero.flows ) {
    print( Amat , quote = FALSE )
    stop( "Some observed flow cells had zero counts.
             If those are zero counts in the population transition probability matrix, consider using as.zero.flows = TRUE.
             If not, consider collapsing categories." )
  } else if ( any( Nij <= 0 ) & as.zero.flows ) {
    warning( "Some observed flow cells had zero counts. Model fitted with zeroes in the population transition probability matrix." , immediate. = TRUE )
  }

  # treat full response
  if ( all( c( Ri , Cj , M ) == 0 ) ) {

    # model fitting
    N <- sum( Nij )
    mfit <-
      list( "observed.counts" = Amat ,
            eta = rowSums( Nij ) / N ,
            pij = sweep( Nij , 1 , rowSums( Nij ) , "/" ) )
    mfit$nipij <- sweep( mfit$pij , 1 , mfit$eta , "*" )
    mfit$gamma <- colSums( mfit$nipij )
    mfit$delta <- N * ( mfit$gamma - mfit$eta )
    mfit$muij <- N*sweep( mfit$pij , 1 , mfit$eta , "*" )

    # visit proportion
    model.expected <- rbind( cbind( mfit$nipij , rep( 0 , ncol( Nij ) ) ) , rep( 0 , nrow( Nij ) + 1 ) )
    interview.number <- data.frame( visit = rep( 1 , length( ww ) ) )
    vVec <- c( stats::xtabs( ww ~ visit , data = interview.number , addNA = TRUE , drop.unused.levels = FALSE ) / sum( ww ) )
    model.expected <- outer( model.expected , vVec )

    # rao-scott adjustment of the chi-square
    K <- dim( model.expected )[1] - 1
    n.parms <- K + K^2
    pearson <- rao.scott( xx , interview.number , ww , model.expected , design , n.parms , pij.zero = sum( mfit$pij == 0 ) )

    # linearization
    llin <- FullResponse.linearization( xx , ww , mfit , design )

    # calculate variance
    mvar <- lapply( llin , function(z) stats::vcov( survey::svytotal( z , design ) ) )

    # build results list
    res <- sapply( c( "eta" , "pij" , "muij" , "gamma" , "delta" ) , function(z) {
      if ( z %in% c( "psi" , "rho" , "tau" , "eta" , "gamma" , "delta" ) ) {
        this.stats <- mfit[[z]]
        attr( this.stats , "var" ) <- mvar[[z]]
        names( this.stats ) <- if ( length( this.stats ) > 1 ) xlevels else z
        colnames( attr( this.stats , "var" ) ) <- rownames( attr( this.stats , "var" ) ) <- if ( length( attr( this.stats , "var" ) ) > 1 ) xlevels else z
        class( this.stats ) <- "svystat"
        attr( this.stats , "statistic" ) <- z
      } else if ( z %in% c( "pij" , "muij" ) ) {
        this.stats <- mfit[[z]]
        these.classes <- expand.grid( dimnames( this.stats ) )
        these.classes <- these.classes[ order( these.classes[, 2 ] ) , ]
        these.classes <- apply( these.classes , 1 , paste, collapse = ":" )
        this.vector <- c( t( this.stats ) )
        names( this.vector ) <- these.classes
        this.vmat <-  mvar[[z]]
        dimnames( this.vmat ) <- list( these.classes , these.classes )
        attr( this.vector , "var" ) <- this.vmat
        attr( this.vector , "categories" ) <- dimnames( this.stats )
        class( this.vector ) <- c( "svymstat" , "svystat" )
        attr( this.vector , "statistic" ) <- z
        this.stats <- this.vector
      }
      return(this.stats)
    } , simplify = FALSE )

    # add influence attribute
    for ( this.stat in names( res ) ) attr( res[[ this.stat ]] , "influence" ) <- if ( influence ) llin[[this.stat]] else NULL

    # create final object
    rval <- res[ c( "eta" , "gamma" , "pij" , "muij" , "delta" ) ]
    rval$model <- "Full Response"
    class(rval) <- "flowstat"
    attr( rval , "formula" ) <- x
    attr( rval , "has.order" ) <- has.order
    attr( rval , "iter" ) <- NA
    attr( rval , "adj.chisq" ) <- pearson

    # return final object
    return( rval )

  }

  # model fitting
  mfit <- ipf( Amat , model , tol = tol , maxit = maxit , verbose = verbose , keep.info = FALSE )

  # visit proportions
  interview.number <- data.frame( visit = rep( 1 , length( ww) ) )
  vVec <- c( stats::xtabs( ww ~ visit , data = interview.number , addNA = TRUE , drop.unused.levels = FALSE ) / sum( ww ) )
  model.expected <- outer( mfit$estimated.props , vVec )

  # rao-scott adjustment of the chi-square
  K <- dim( model.expected )[1] - 1
  n.parms <- switch (model, A = { K^2 + K + 3 } , B = { K^2 + 2*K + 2 } , C = { K^2 + 3*K + 1 } , D = { K^2 + 3*K + 1 } )
  pearson <- rao.scott( xx , interview.number , ww , model.expected , design , n.parms , pij.zero = sum( mfit$pij == 0 ) )

  # estimate linearized variables
  llin <- linearization_fun( xx , ww , res = mfit , design = design )
  llin[ lapply( llin , class ) %in% "numeric" ] <- lapply( llin[ lapply( llin , class ) %in% "numeric" ] , matrix , nrow = length( ww ) )

  # calculate variance
  mvar <- lapply( llin , function(z) stats::vcov( survey::svytotal( z , design ) ) )

  # build results list
  res <- sapply( c( "psi" , "rho" , "tau" , "eta" , "pij" , "muij" , "gamma" , "delta" ) , function(z) {
    if ( z %in% c( "psi" , "rho" , "tau" , "eta" , "gamma" , "delta" ) ) {
      this.stats <- mfit[[z]]
      attr( this.stats , "var" ) <- mvar[[z]]
      names( this.stats ) <- if ( length( this.stats ) > 1 ) xlevels else z
      colnames( attr( this.stats , "var" ) ) <- rownames( attr( this.stats , "var" ) ) <- if ( length( attr( this.stats , "var" ) ) > 1 ) xlevels else z
      class( this.stats ) <- "svystat"
      attr( this.stats , "statistic" ) <- z
    } else if ( z %in% c( "pij" , "muij" ) ) {
      this.stats <- mfit[[z]]
      these.classes <- expand.grid( dimnames( this.stats ) )
      these.classes <- these.classes[ order( these.classes[, 2 ] ) , ]
      these.classes <- apply( these.classes , 1 , paste, collapse = ":" )
      this.vector <- c( t( this.stats ) )
      names( this.vector ) <- these.classes
      this.vmat <-  mvar[[z]]
      dimnames( this.vmat ) <- list( these.classes , these.classes )
      attr( this.vector , "var" ) <- this.vmat
      attr( this.vector , "categories" ) <- dimnames( this.stats )
      class( this.vector ) <- c( "svymstat" , "svystat" )
      attr( this.vector , "statistic" ) <- z
      this.stats <- this.vector
    }
    return(this.stats)
  } , simplify = FALSE )

  # add influence attribute
  for ( this.stat in names( res ) ) attr( res[[ this.stat ]] , "influence" ) <- if ( influence ) llin[[this.stat]] else NULL

  # create final object
  rval <- res[ c( "psi" , "rho" , "tau" , "eta" , "gamma" , "pij" , "muij" , "delta" ) ]
  rval$model <- mfit$model
  class(rval) <- "flowstat"
  attr( rval , "formula" ) <- x
  attr( rval , "has.order" ) <- has.order
  attr( rval , "iter" ) <- mfit$iter
  attr( rval , "adj.chisq" ) <- pearson

  # return final object
  rval

}

#' @export
svyflow <- function( x , design , model = c( "A","B","C","D") , ... ) {

  # test values
  model <- match.arg( model , several.ok = FALSE )

  # test valid arguments
  if ( class( x ) != "formula" ) stop( "x must be a formula." )
  if ( length( as.character( x ) ) != 2 ) stop( "x must be a one-sided formula." )
  if ( length( ( strsplit( as.character( x )[[2]] , " \\+ " ) )[[1]] ) != 2 ) stop("only two-way tables at the moment.")
  if ( ncol( attr( terms( x ) , "factors" ) ) != 2 ) stop("only two-way tables at the moment.")
  UseMethod( "svyflow" , design )

}
