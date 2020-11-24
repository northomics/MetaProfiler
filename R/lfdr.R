#' Computes the default argument of various distributions. Adapted from fitdistrplus.
#'
#' @param x A numeric vector
#' @param distr The name of the distribution. Can be either "weibull", "beta", "norm", 
#' "lnorm", "exp", "gamma", "invtrgamma", "trgamma", "lgamma", "pareto", "pareto1", 
#' "invweibull", "llogis", "invgamma", "unif", "cauchy", or "logis".
#' @return A list of starting values for \code{distr}
#' @examples
#' x <- rnorm(100)
#' args <- start_arg_default(x, "norm")
#' args
start_arg_default <- function (x, distr)
{
  if (distr == "norm") {
    n <- length(x)
    sd0 <- sqrt((n - 1)/n) * sd(x)
    mx <- mean(x)
    start <- list(mean = mx, sd = sd0)
  }
  else if (distr == "lnorm") {
    if (any(x <= 0)) 
      stop("values must be positive to fit a lognormal distribution")
    n <- length(x)
    lx <- log(x)
    sd0 <- sqrt((n - 1)/n) * sd(lx)
    ml <- mean(lx)
    start <- list(meanlog = ml, sdlog = sd0)
  }
  else if (distr == "exp") {
    if (any(x < 0)) 
      stop("values must be positive to fit an exponential  distribution")
    start <- list(rate = 1/mean(x))
  }
  else if (distr == "gamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an gamma  distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n * var(x)
    start <- list(shape = m^2/v, rate = m/v)
  }
  else if (distr == "beta") {
    if (any(x < 0) | any(x > 1)) 
      stop("values must be in [0-1] to fit a beta distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n * var(x)
    aux <- m * (1 - m)/v - 1
    start <- list(shape1 = m * aux, shape2 = (1 - m) * aux)
  }
  else if (distr == "weibull") {
    if (any(x < 0)) 
      stop("values must be positive to fit an Weibull  distribution")
    m <- mean(log(x))
    v <- var(log(x))
    shape <- 1.2/sqrt(v)
    scale <- exp(m + 0.572/shape)
    start <- list(shape = shape, scale = scale)
  }
  else if (distr == "logis") {
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n * var(x)
    start <- list(location = m, scale = sqrt(3 * v)/pi)
  }
  else if (distr == "cauchy") {
    start <- list(location = median(x), scale = IQR(x)/2)
  }
  else if (distr == "unif") {
    start <- list(min = 0, max = 1)
  }
  else if (distr == "invgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse gamma  distribution")
    m1 <- mean(x)
    m2 <- mean(x^2)
    shape <- (2 * m2 - m1^2)/(m2 - m1^2)
    scale <- m1 * m2/(m2 - m1^2)
    start <- list(shape = shape, scale = scale)
  }
  else if (distr == "llogis") {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-logistic  distribution")
    p25 <- as.numeric(quantile(x, 0.25))
    p75 <- as.numeric(quantile(x, 0.75))
    shape <- 2 * log(3)/(log(p75) - log(p25))
    scale <- exp(log(p75) + log(p25))/2
    start <- list(shape = shape, scale = scale)
  }
  else if (distr == "invweibull") {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse Weibull distribution")
    g <- log(log(4))/(log(log(4/3)))
    p25 <- as.numeric(quantile(x, 0.25))
    p75 <- as.numeric(quantile(x, 0.75))
    shape <- exp((g * log(p75) - log(p25))/(g - 1))
    scale <- log(log(4))/(log(shape) - log(p25))
    start <- list(shape = shape, scale = max(scale, 1e-09))
  }
  else if (distr == "pareto1") {
    if (any(x < 0)) 
      stop("values must be positive to fit a Pareto distribution")
    x1 <- min(x)
    m1 <- mean(x)
    n <- length(x)
    shape <- (n * m1 - x1)/(n * (m1 - x1))
    min <- x1 * (n * shape - 1)/(n * shape)
    start <- list(shape = shape, min = min)
  }
  else if (distr == "pareto") {
    if (any(x < 0)) 
      stop("values must be positive to fit a Pareto distribution")
    m1 <- mean(x)
    m2 <- mean(x^2)
    scale <- (m1 * m2)/(m2 - 2 * m1^2)
    shape <- 2 * (m2 - m1^2)/(m2 - 2 * m1^2)
    start <- list(shape = shape, scale = scale)
  }
  else if (distr == "lgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-gamma distribution")
    m1 <- mean(log(x))
    m2 <- mean(log(x)^2)
    alpha <- m1^2/(m2 - m1^2)
    lambda <- m1/(m2 - m1^2)
    start <- list(shapelog = alpha, ratelog = lambda)
  }
  else if (distr == "trgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an trans-gamma  distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n * var(x)
    start <- list(shape1 = m^2/v, shape2 = 1, rate = m/v)
  }
  else if (distr == "invtrgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse trans-gamma  distribution")
    n <- length(1/x)
    m <- mean(1/x)
    v <- (n - 1)/n * var(1/x)
    start <- list(shape1 = m^2/v, shape2 = 1, rate = m/v)
  }
  else stop(paste0("Unknown starting values for distribution ", distr, "."))
  return(start)
}

.fit_true_distribution <- function(data, n_iter, kde0, distr, times, start_prior_prob, trace, ...) {
  if(trace) {
    on.exit({
      sink()
      close(con)
    })
  }
  start_arg = mapply(function(dname, dat) {
    start_arg = start_arg_default(dat, dname)
  }, distr, data, SIMPLIFY = FALSE)
  data1 = data
  if(n_iter < nrow(data)) {
    i = sample(seq_len(nrow(data)), n_iter, replace = TRUE)
    data1 = data[i,]
  }
  u1 = VineCopula::pobs(data1)
  cop = copula::normalCopula(dim = ncol(data))
  fit <- copula::fitCopula(cop, u1)
  estimates = fit@estimate
  vstart = c(estimates, unlist(start_arg), start_prior_prob)
  nestim = length(estimates)
  split_by = rep(seq_len(ncol(data)), lengths(start_arg))
  arg_names = unlist(lapply(start_arg, names))
  obs = as.matrix(data1)
  nparams = length(vstart) - 1
  p0 = ks::dkde(obs, kde0)
  stdout <- character()
  envir <- environment()
  start <- 1
  mle = function(x) {
    args = setNames(x[(nestim+1):nparams], arg_names)
    args = split(setNames(as.list(args), arg_names), split_by)
    cop@parameters = x[1:nestim]
    mv_cop = try(suppressMessages(copula::mvdc(cop, distr, args)), silent = TRUE)
    if(inherits(mv_cop, "try-error")) {
      return(.Machine$double.xmax)
    }
    prior_prob = x[nparams + 1]
    if (prior_prob < 0 || prior_prob > 1) {
      return(.Machine$double.xmax)
    }
    p1 = try(copula::dMvdc(obs, mv_cop), silent = TRUE)
    if(inherits(p1, "try-error")) {
      return(.Machine$double.xmax)
    }
    p1[is.nan(p1)] = 0
    p = p0 * prior_prob + p1 * (1 - prior_prob)
    val = -sum(log(p))
    if(trace) {
      end = length(stdout)
      if(start <= end) {
        message(crayon::silver(paste0(stdout[start:end], collapse = "\n")))
      }
      assign("start", length(stdout) + 1, envir = envir)
    }
    val
  }
  a <- list(...)
  if(!any(names(a) == "method")) {
    a <- c(a, list(method = "Nelder-Mead"))
  }
  if(trace) {
    con <- textConnection('stdout', 'wr', local = TRUE)
    sink(con)
  }
  opt <- do.call(optim, c(list(par = vstart, f = mle), a))
  if(length(a$control$trace) && a$control$trace) {
    end = length(stdout)
    message(crayon::silver(paste0(stdout[start:end], collapse = "\n")))
  }
  x = opt$par
  prior_prob = x[nparams + 1]
  args = setNames(x[(nestim+1):nparams], arg_names)
  args = split(setNames(as.list(args), arg_names), split_by)
  cop@parameters = x[1:nestim]
  mv_cop = suppressMessages(copula::mvdc(cop, distr, args))
  list(copula = cop, mvdc = mv_cop, prior_prob = prior_prob)
}

#' Computes the LFDR of heavy features quantified from Protein-SIP experiments.
#'
#' @param Object An object of class MetaProfiler.
#' @param condition_names A character vector specifying the conditions used in the experiment. 
#' Names must match the columns in Object@@design. If \code{NULL}, then all the conditions are 
#' combined to compute LFDR.
#' @param estimate_prior_probability Either "mle" or "analytically". If "mle" is chosen, prior 
#' probability for false discovery, p0, is estimated from \eqn{f(x) = f0(x)*p0 + f1(x)*(p0-1)} 
#' using maximum likelihood estimation, where \eqn{f0(x)}, \eqn{f1(x)}, and \eqn{f(x)} is the 
#' density of the data for false, true, and mixed discoveries, respectively. If "analytically" 
#' is chosen, prior probability is simply the number of features at time 0 of the corresponding 
#' \code{by} group divided by the number of features in the current group.
#' @param distr (Only used when \code{estimate_prior_probability} is set to \code{"mle"}) 
#' A vector with two strings values. The first specifies the name of the distribution of 
#' RIA for true discoveries and the second specifies the distribution of LR for true discoveries. 
#' Can be either "weibull", "beta", "norm", "lnorm", "exp", "gamma", "invtrgamma", "trgamma", 
#' "lgamma", "pareto", "pareto1", "invweibull", "llogis", "invgamma", "unif", "cauchy", or "logis".
#' @param n_iter (Only used when \code{estimate_prior_probability} is set to \code{"mle"}) 
#' A numeric value for the number of iterations used for the bootstrap. The bootstrap values are 
#' used to calculate the log likelihood of the parameters for the normal copula.
#' @param trace A logical value. If \code{TRUE}, then log is printed. 
#' @param progress A function. One over the total number of loops that will be performed
#' will be passed as the first argument. The total number of loops should be equal to the total
#' number of timepoints times, if \code{condition_names} is not \code{NULL}, 
#' the number of conditions.
#' @param ... Additional arguments to be passed to \code{progress}.
#' @return Returns an object of class MetaProfiler.
#' @examples
lfdr <- function(Object, observations = c(Object@incorporation_name, Object@labeling_ratio_name),
                 condition_names = character(0), estimate_prior_probability = c("mle", "analytically"),
                 distr = rep("weibull", 2), n_iter = 10000, trace = TRUE, progress = NULL,
                 score_threshold = NULL, id_type = c("both", "id", "feature"), seed = 123, ...) {
  set.seed(seed)
  if(length(Object@data$Feature)) {
    id_type = match.arg(id_type)
    if(id_type != "both") {
      Object@data[Feature == id_type]
    }
  }
  if (is.function(progress) && 
      !all(deparse(progress) == deparse(shiny::incProgress))) {
    trace <- FALSE
  }
  estimate_prior_probability = match.arg(estimate_prior_probability, c("mle", "analytically"))
  dname = distr
  ddname = paste0("d", dname)
  viable_ddname = exists(ddname, mode = "function")
  if(!viable_ddname) {
    n = sum(!viable_ddname)
    msg = combine_words(ddname[!viable_ddname])
    stop("The function ", ddname, " is not defined.")
  }
  data = Object@data
  if(length(score_threshold)) {
    if(Object@higher_score_better)
      Object@data = data[get(Object@score_name) >= score_threshold]
    else
      Object@data = data[get(Object@score_name) <= score_threshold]
  }
  if (length(condition_names))
    Object@data[,condition_names] <- Object@data[,lapply(.SD, as.character), .SDcols = condition_names]
  by = c(condition_names, Object@time_unit)
  design <- unique(Object@design[get(Object@time_unit) != Object@time_zero, ..by])
  design <- design[order(get(Object@time_unit))]
  
  
  design0 <- unique(Object@design[get(Object@time_unit) == Object@time_zero, ..by])
  design0 <- design0[order(get(Object@time_unit))]
  
  if(!length(Object@data$N)) {
    Object@data$N = 1
  }
  false <- Object@data[get(Object@time_unit) == Object@time_zero]
  false_observation_table <- false[,c(condition_names, observations), with = F]
  kde0s = false_observation_table[, .(pdf = list(
    list(copula = NULL, mvdc = NULL, prior_prob = 1, data = .SD, kde = NULL, kde0 = ks::kde(as.matrix(.SD))))
  ), by = condition_names]
  if(length(condition_names)) {
    kde0s = kde0s[design0[, ..condition_names], ,  on = condition_names]
  }
  Object@pdf[apply(design0, 1, function(x) paste(by, x, collapse = ", "))] = kde0s$pdf
  Object@data$id <- seq_len(nrow(Object@data))
  Object@data$LFDR <- NA_real_
  for(i in 1:nrow(design)) {
    mixture <- Object@data[design[i,], , on = by]
    mixture_observation_table <- mixture[, ..observations]
    if(length(condition_names)){
      kde0 = kde0s[unique(mixture[,..condition_names]), , on = condition_names]$pdf[[1]]$kde0
    } else {
      kde0 = kde0s$pdf[[1]]$kde0
    }
    prior_prob = sum(false$N)/sum(mixture$N)
    if (prior_prob > 1) {
      prior_prob <- 0.8
    }
    x = as.matrix(mixture_observation_table[, lapply(.SD, function(x) rep(x, mixture$N))])
    kde = ks::kde(x)
    if(estimate_prior_probability == "mle" & length(distr)) {
      if (trace) {
        cat(crayon::blue("Fitting ", dplyr::case_when(length(observations) == 1 ~ "uni",
                                                      length(observations) == 2 ~ "bi",
                                                      length(observations) == 3 ~ "tri",
                                                      T ~ "multi"),"variate Gaussian copula to ",
                         combine_words(observations),
                         " for true discoveries at ", tolower(Object@time_unit), " ",
                         unique(mixture[,get(Object@time_unit)]), 
                         ifelse(length(condition_names), list(paste0(
                           " in ", design[i,..condition_names])), list(NULL))[[1]], "...\n", sep = ""))
        cat(crayon::silver("Initial prior probability for false discoveries: ", round(prior_prob * 100), "%\n", sep = ""))
      }
      fit <- .fit_true_distribution(mixture_observation_table, n_iter, kde0, distr, mixture$N, prior_prob, trace, ...)
      x = as.matrix(mixture_observation_table)
      prior_prob = fit$prior_prob
      if (trace) {
        cat(crayon::silver("Final prior probability for false discoveries:", round(prior_prob * 100), "%\n", sep = ""))
      }
      odds = (prior_prob)/(1 - prior_prob) * (ks::dkde(x, kde0)/copula::dMvdc(x, fit$mvdc))
      LFDR = 1/(1 + odds^-1)
    } else {
      fit = NULL
      LFDR = prior_prob * ks::dkde(x, kde0)/ks::dkde(x, kde)
    }
    Object@data$LFDR[mixture$id] = LFDR
    Object@pdf[paste(by, design[i,], collapse = ", ")] = list(c(fit, list(data = mixture_observation_table, kde = kde, kde0 = kde0)))
    if(is.function(progress)) {
      progress(1/nrow(design))
    }
  }
  Object
}


#' Filter data according to thresholds.
#'
#' @param Object An object of class MetaProfiler.
#' @param LFDR_threshold Either a number or a character value "strong" or "substantial". 
#' If "strong" is chosen, then the threshold is 0.0909 (odds for false discoveries is 1 out of 10). 
#' If "substantial" is chosen, then the LFDR threshold is 0.2 (odds for false discovery is 1 out 
#' of 3). 
#' @param min_nb_timepoints Numeric value specifying the minimum number of timepoints the peptide or 
#' protein must be present in.
#' @param min_labeling_ratio (Only used when \code{estimate_prior_probability} is set to "mle") A vector 
#' with two strings values. The first specifies the name of the distribution of RIA for true discoveries 
#' and the second specifies the distribution of LR for true discoveries. Can be either "weibull", 
#' "beta", "norm", "lnorm", "exp", "gamma", "invtrgamma", "trgamma", "lgamma", "pareto", 
#' "pareto1", "invweibull", "llogis", "invgamma", "unif", "cauchy", or "logis".
#' @param condition_names A character vector specifying the conditions used in the experiment. 
#' Names must match the columns in Object@@design. Defaults to all columns in Object@@design, 
#' except for the column containing the timepoints. Names must match the columns in Object@@design.
#' @param min_nb_samples Numeric value specifying the minimum number of conditions the peptide or protein 
#' must be present in.
#' @param score_threshold The score threshold.
#' @param higher_score_better Logical value denoting if higher score is better.
#' @param trace If \code{TRUE}, then progress is printed. Defaults to \code{TRUE}.
#' @return Returns an object of class MetaProfiler.
#' @examples
#' 

filter_data <- function(Object, 
                        LFDR_threshold = c("strong", "substantial"),
                        score_threshold = NULL,
                        min_nb_timepoints = NULL, 
                        condition_names = colnames(Object@design)[colnames(Object@design) != Object@time_unit],
                        min_nb_samples = NULL,
                        combine_samples = T,
                        aggregate = c("weighted", "shrink", "mean", "median"),
                        id_type = c("both", "feature", "id"),
                        clust_features = T,
                        timepoints = Object@timepoints,
                        ...,
                        trace = TRUE) {
  id_type = match.arg(id_type)
  if(clust_features) {
    data <- cluster_features(Object, id_type = id_type, ...)@data[get(Object@time_unit) %in% timepoints]
  } else {
    if(id_type != "both" && length(Object@data$Feature)) {
      data <- data.table::copy(Object@data[get(Object@time_unit) %in% timepoints & Feature == id_type])
    } else {
      data <- data.table::copy(Object@data[get(Object@time_unit) %in% timepoints])
    }
  }
  by = c(Object@peptide_column_PTMs, Object@peptide_column_no_PTMs, Object@accession_column)
  if(!combine_samples && length(condition_names)) {
    by = c(by, condition_names)
  }
  data = data[get(Object@time_unit) != Object@time_zero]
  if(LFDR_threshold[1] == "substantial") {
    LFDR_threshold <- 1-1/(1+1/3)
  } else if (LFDR_threshold[1] == "strong") {
    LFDR_threshold <- 1-1/(1+1/10)
  }
  data <- data[LFDR <= LFDR_threshold]
  if(length(min_nb_samples) && min_nb_samples > 0 && length(condition_names)) {
    data <- data[,test := (length(unique(get(Object@time_unit))) >= min_nb_timepoints) & (nrow(unique(.SD)) >= min_nb_samples), by = by, .SDcols = condition_names]
  } else {
    min_nb_timepoints <- ifelse(length(min_nb_timepoints), min_nb_timepoints, 0)
    data <- data[,test := (length(unique(get(Object@time_unit))) >= min_nb_timepoints), by = by]
  }
  data <- data[(test),-"test"]
  if(trace) {
    message('***** LFDR Filtering *****')
    cat(length(unique(data[[Object@peptide_column_PTMs]])),
        " peptides left after filtering at an LFDR of ",
        LFDR_threshold * 100, "%.\n", sep = "")
  }
  
  col <- na.omit(c(Object@incorporation_columns[1], Object@incorporation_name[1], Object@intensity_columns[1], Object@intensity_name[1], Object@score_columns[1], Object@score_name[1], Object@labeling_ratio_name[1], "LFDR"))
  LR <- data[,get(Object@labeling_ratio_name)]
  if(.check_if_percantage_or_fraction(LR)) {
    LR <- LR/100
  }
  if(!is.na(Object@intensity_columns[1]))
    data[,Object@intensity_name] <- LR * data[[Object@intensity_columns[1]]]/(1 - LR)
  dt <- .aggregate_by(data = data, by = by, tu = Object@time_unit, 
                      vars = c(Object@incorporation_name[1], Object@labeling_ratio_name[1]),
                      col = col, method = aggregate, trace = trace)

  by2 = by
  by2[!(by2 %in% make.names(by2))] <- paste0('`', by2[!(by2 %in% make.names(by2))], '`')
  formula = eval(parse(text = paste(paste0(by2, collapse = "+"), " ~ ", Object@time_unit)))

  dt <- data.table::dcast(dt, formula = formula, value.var = col)
  timepoints <- c(Object@time_zero, timepoints)
  col_add <- c(by, paste0(rep(col, each = length(timepoints)), "_", rep(timepoints, length(col))))
  dt[,setdiff(col_add, colnames(dt))] <- NA
  dt <- dt[,..col_add]
  new_cols = c(paste("Light", Object@incorporation_name),
               paste("Heavy", Object@incorporation_name),
               paste("Light", Object@intensity_name),
               paste("Heavy", Object@intensity_name),
               paste("Light", Object@score_name),
               paste("Heavy", Object@score_name),
               Object@labeling_ratio_name, "LFDR")
  if(length(attr(col, "na.action")))
    new_cols <- new_cols[-(attr(col, "na.action"))]
  colnames(dt)[-(seq_along(by))] <- get_cols(Object, new_cols, timepoints)
  Object@master_tbl <- dt
  if(combine_samples)
    condition_names = NULL
  Object <- get_annotations(Object, condition_names)
  Object
}
