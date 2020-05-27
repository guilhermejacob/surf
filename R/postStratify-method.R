#' @method postStratify surflow.design
#' @export
postStratify.surflow.design <-
  function (design, strata, population, partial = FALSE, ...) {
  if (inherits(strata, "formula")) {
    mf <- substitute(model.frame(strata, data = design$variables[[1]], na.action = na.fail))
    strata <- eval.parent(mf)
  }
  strata <- as.data.frame(strata)
  sampletable <- stats::xtabs(I(1/design$prob) ~ ., data = strata)
  sampletable <- as.data.frame(sampletable)
  if (inherits(population, "table")) population <- as.data.frame(population)
  else if (is.data.frame(population)) population$Freq <- as.vector(population$Freq) else stop("population must be a table or dataframe")
  if (!all(names(strata) %in% names(population)))
    stop("Stratifying variables don't match")
  nn <- names(population) %in% names(strata)
  if (sum(!nn) != 1)
    stop("stratifying variables don't match")
  names(population)[which(!nn)] <- "Pop.Freq"
  both <- merge(sampletable, population, by = names(strata),
                all = TRUE)
  samplezero <- both$Freq %in% c(0, NA)
  popzero <- both$Pop.Freq %in% c(0, NA)
  both <- both[!(samplezero & popzero), ]
  if (any(onlysample <- popzero & !samplezero)) {
    print(both[onlysample, ])
    stop("Strata in sample absent from population. This Can't Happen")
  }
  if (any(onlypop <- samplezero & !popzero)) {
    if (partial) {
      both <- both[!onlypop, ]
      warning("Some strata absent from sample: ignored")
    }
    else {
      print(both[onlypop, ])
      stop("Some strata absent from sample: use partial=TRUE to ignore them.")
    }
  }
  reweight <- both$Pop.Freq/both$Freq
  both$label <- do.call("interaction", list(both[, names(strata)]))
  designlabel <- do.call("interaction", strata)
  index <- match(designlabel, both$label)
  attr(index, "oldweights") <- 1/design$prob
  design$prob <- design$prob/reweight[index]
  attr(index, "weights") <- 1/design$prob
  design$postStrata <- c(design$postStrata, list(index))
  design$call <- sys.call(-1)
  design
}
