DEoptim.control <- function (VTR = -Inf, strategy = 2, bs = FALSE, NP = 50, itermax = 200, 
                             CR = 0.5, F = 0.8, trace = TRUE, initialpop = NULL, storepopfrom = itermax + 1,
                             storepopfreq = 1, p = 0.2, c = 0, reltol = sqrt(.Machine$double.eps), 
                             steptol = itermax) 
{
  if (itermax <= 0) {
    warning("'itermax' <= 0; set to default value 200\n", 
            immediate. = TRUE)
    itermax <- 200
  }
  if (NP < 4) {
    warning("'NP' < 4; set to default value 50\n", immediate. = TRUE)
    NP <- 50
  }
  if (F < 0 | F > 2) {
    warning("'F' not in [0,2]; set to default value 0.8\n", 
            immediate. = TRUE)
    F <- 0.8
  }
  if (CR < 0 | CR > 1) {
    warning("'CR' not in [0,1]; set to default value 0.5\n", 
            immediate. = TRUE)
    CR <- 0.5
  }
  if (strategy < 1 | strategy > 6) {
    warning("'strategy' not in {1,...,6}; set to default value 2\n", 
            immediate. = TRUE)
    strategy <- 2
  }
  bs <- (bs > 0)
  if (trace < 0) {
    warning("'trace' cannot be negative; set to 'TRUE'")
    trace <- TRUE
  }
  storepopfreq <- floor(storepopfreq)
  if (storepopfreq > itermax) 
    storepopfreq <- 1
  if (p <= 0 || p > 1) {
    warning("'p' not in (0,1]; set to default value 0.2\n", 
            immediate. = TRUE)
    p <- 0.2
  }
  if (c < 0 || c > 1) {
    warning("'c' not in [0,1]; set to default value 0\n", 
            immediate. = TRUE)
    c <- 0
  }
  list(VTR = VTR, strategy = strategy, NP = NP, itermax = itermax, 
       CR = CR, F = F, bs = bs, trace = trace, initialpop = initialpop, 
       storepopfrom = storepopfrom, storepopfreq = storepopfreq, 
       p = p, c = c, reltol = reltol, steptol = steptol)
}

require(Rcpp)

sourceCpp("~MetaProfiler/curve_fitting.cpp")

curve_fitting <- function(Object, var, use_knn = T, lower = c(0,0,0,0),
                          upper = c(kbi = 100,100,100,100), control = DEoptim.control(trace = F),
                          equation = 3, progress = T, seed = 362436069,
                          timepoints = Object@timepoints,
                          data = NULL) {
  return_consts_only = T
  if(!length(data)) {
    data = Object@master_tbl
    return_consts_only = F
  }
  # if(is.null(fn)) {
  #   fn = cppFunction(
  #     'double three_exponential_function(const Rcpp::NumericVector & x, const Rcpp::NumericVector & xr, const Rcpp::NumericVector & yr) {
  #       double u = ((x[3] + x[1] + x[3]) - std::sqrt(std::pow((x[3] + x[1] + x[3]), 2) - (4 * x[1]*x[3]))) / 2;
  #       double v = ((x[3] + x[1] + x[3]) + std::sqrt(std::pow((x[3] + x[1] + x[3]), 2) - (4 * x[1]*x[3]))) / 2;
  #       double yu = x[1] * x[2]*(u - x[3]) / ((u - v)*(u - x[2])*u);
  #       double yv = x[1] * x[2]*(v - x[3]) / ((v - u)*(v - x[2])*v);
  #       double ykbi = x[1] * (x[2] - x[3]) / ((u - x[2])*(v - x[2]));
  #       double sum = 0;
  #       for(int i; i < xr.size(); ++i) {
  #         double value = (1 + yu * exp(-u * xr[i]) + yv * exp(-v * xr[i]) + ykbi * exp(-x[2] * xr[i]))*100;
  #         sum += std::pow(value - yr[i], 2);
  #       }
  #       return sum;
  #     }'
  #   )
  # }
  ##fn1  <- function(par) fn(par, ...)
  if (length(lower) != length(upper))
    stop("'lower' and 'upper' are not of same length")
  if (!is.vector(lower))
    lower <- as.vector(lower)
  if (!is.vector(upper))
    upper <- as.vector(upper)
  if (any(lower > upper))
    stop("'lower' > 'upper'")
  if (any(lower == "Inf"))
    warning("you set a component of 'lower' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(lower == "-Inf"))
    warning("you set a component of 'lower' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "Inf"))
    warning("you set a component of 'upper' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "-Inf"))
    warning("you set a component of 'upper' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (!is.null(names(lower))) {
    consts <- names(lower)
    params <- consts
  }
  else if (!is.null(names(upper)) & is.null(names(lower))){
    consts <- names(upper)
    params <- consts
  }
  else {
    consts <- paste("par", 1:length(lower), sep = "")
    params <- consts
  }
  consts <- c(consts, "score")
  col <- get_cols(Object, var, timepoints)
  mat <- data[,..col,with=F]
  if(use_knn)
  {
    shush(mat <- data.table::as.data.table(
      impute::impute.knn(as.matrix(mat), rowmax = 1, colmax = 1)$data
    ))
  }
  
  mat <- data.table::data.table(
    xr = list(timepoints),
    yr = unlist(apply(as.matrix(mat), 1, function(x) list(x)), recursive = F)
  )
  # makeEnv <- function(...) {
  #   list2env(list(...))
  # }
  # env <- mapply(function(xr, yr){
  #   makeEnv(xr = xr, yr = yr)
  # }, mat$xr, mat$yr)
  ctrl <- do.call(DEoptim.control, as.list(control))
  ctrl$npar <- length(lower)
  if (ctrl$NP < 4) {
    warning("'NP' < 4; set to default value 50\n", immediate. = TRUE)
    ctrl$NP <- 50
  }
  if (ctrl$NP < 10*length(lower))
    warning("For many problems it is best to set 'NP' (in 'control') to be at least ten times the length of the parameter vector. \n", immediate. = TRUE)
  if (!is.null(ctrl$initialpop)) {
    ctrl$specinitialpop <- TRUE
    if(!identical(as.numeric(dim(ctrl$initialpop)), c(ctrl$NP, ctrl$npar)))
      stop("Initial population is not a matrix with dim. NP x length(upper).")
  } else {
    ctrl$specinitialpop <- FALSE
    ctrl$initialpop <- matrix(0,1,1)    # dummy matrix
  }
  ##
  ctrl$trace <- as.numeric(ctrl$trace)
  ctrl$specinitialpop <- as.numeric(ctrl$specinitialpop)
  set.seed(seed)
  mat[, consts] <- data.table::as.data.table(t(curve_fitting_c(mat$xr,
                                                  mat$yr,
                                                  minbound = lower,
                                                  maxbound = upper,
                                                  equation = equation,
                                                  control = ctrl,
                                                  verbose = progress)))
  mat$equation <- equation
  consts_tbl <- mat[,c(consts, "equation"), with = F]
  consts_tbl[,paste0(params, "_lower")] <- data.table::as.data.table(
    matrix(lower, nrow = nrow(consts_tbl),
           ncol = length(params), byrow = T)
  )
  consts_tbl[,paste0(params, "_upper")] <- data.table::as.data.table(
    matrix(upper, nrow = nrow(consts_tbl),
           ncol = length(params), byrow = T)
  )
  if(return_consts_only) {
    return(consts_tbl)
  }
  Object@consts[var] <- list(consts_tbl)
  Object
}

validate_curve_fitting <- function(Object, var,
                               use_knn = T,
                               lower = c(0,0,0,0),
                               upper = c(100,100,100,100),
                               control = DEoptim.control(trace = F),
                               progress = T,
                               equation = 3,
                               timepoints = Object@timepoints,
                               seed = 362436069, data = NULL) {
  if(!length(data)) {
    data = Object@master_tbl
  }
  # if(is.null(fn)) {
  #   fn = cppFunction(
  #     'double three_exponential_function(const Rcpp::NumericVector & x, const Rcpp::NumericVector & xr, const Rcpp::NumericVector & yr) {
  #       double u = ((x[3] + x[1] + x[3]) - std::sqrt(std::pow((x[3] + x[1] + x[3]), 2) - (4 * x[1]*x[3]))) / 2;
  #       double v = ((x[3] + x[1] + x[3]) + std::sqrt(std::pow((x[3] + x[1] + x[3]), 2) - (4 * x[1]*x[3]))) / 2;
  #       double yu = x[1] * x[2]*(u - x[3]) / ((u - v)*(u - x[2])*u);
  #       double yv = x[1] * x[2]*(v - x[3]) / ((v - u)*(v - x[2])*v);
  #       double ykbi = x[1] * (x[2] - x[3]) / ((u - x[2])*(v - x[2]));
  #       double sum = 0;
  #       for(int i; i < xr.size(); ++i) {
  #         double value = (1 + yu * exp(-u * xr[i]) + yv * exp(-v * xr[i]) + ykbi * exp(-x[2] * xr[i]))*100;
  #         sum += std::pow(value - yr[i], 2);
  #       }
  #       return sum;
  #     }'
  #   )
  # }
  ##fn1  <- function(par) fn(par, ...)
  if (length(lower) != length(upper))
    stop("'lower' and 'upper' are not of same length")
  if (!is.vector(lower))
    lower <- as.vector(lower)
  if (!is.vector(upper))
    upper <- as.vector(upper)
  if (any(lower > upper))
    stop("'lower' > 'upper'")
  if (any(lower == "Inf"))
    warning("you set a component of 'lower' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(lower == "-Inf"))
    warning("you set a component of 'lower' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "Inf"))
    warning("you set a component of 'upper' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "-Inf"))
    warning("you set a component of 'upper' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  group <- colnames(data)[1]
  col <- get_cols(Object, var, timepoints)
  mat <- data[,c(group, col),with=F]
  if(use_knn) {
    shush(mat[,colnames(mat)[-seq_along(group)]] <- data.table::as.data.table(
      impute::impute.knn(as.matrix(mat[,-..group]), rowmax = 1, colmax = 1)$data
    ))
  }
  mat <- cbind(mat[,..group], data.table::data.table(
    xr = list(timepoints),
    yr = unlist(apply(as.matrix(mat[,-..group]), 1, function(x) list(x)), recursive = F)
  ))
  # makeEnv <- function(...) {
  #   list2env(list(...))
  # }
  # env <- mapply(function(xr, yr){
  #   makeEnv(xr = xr, yr = yr)
  # }, mat$xr, mat$yr)
  ctrl <- do.call(DEoptim.control, as.list(control))
  ctrl$npar <- length(lower)
  if (ctrl$NP < 4) {
    warning("'NP' < 4; set to default value 50\n", immediate. = TRUE)
    ctrl$NP <- 50
  }
  if (ctrl$NP < 10*length(lower))
    warning("For many problems it is best to set 'NP' (in 'control') to be at least ten times the length of the parameter vector. \n", immediate. = TRUE)
  if (!is.null(ctrl$initialpop)) {
    ctrl$specinitialpop <- TRUE
    if(!identical(as.numeric(dim(ctrl$initialpop)), c(ctrl$NP, ctrl$npar)))
      stop("Initial population is not a matrix with dim. NP x length(upper).")
  } else {
    ctrl$specinitialpop <- FALSE
    ctrl$initialpop <- matrix(0,1,1)    # dummy matrix
  }
  ##
  ctrl$trace <- as.numeric(ctrl$trace)
  ctrl$specinitialpop <- as.numeric(ctrl$specinitialpop)
  set.seed(seed)
  ctrl$timepoints <- timepoints
  SE <-curve_fitting_test_c(mat$xr,
                            mat$yr,
                            minbound = lower,
                            maxbound = upper,
                            equation = equation,
                            control = ctrl,
                            verbose = progress)
  SE[t(is.na(Object@master_tbl[,..col]))] <- NA
  RMSE <- sqrt(colMeans(SE, na.rm = T))
}


impute <- function(Object, vars, timepoints = Object@timepoints)
{
  for(var in vars) {
    value = model(Object, var = var, x = timepoints)
    ori = Object@master_tbl[, get_cols(Object, var, timepoints), with = F]
    ori[is.na(ori)] <- value[is.na(ori)]
    Object@master_tbl[,get_cols(Object, var, timepoints)] <- ori
  }
  Object
}


impute2 <- function(Object, vars, timepoints = Object@timepoints)
{
  for(var in impute) {
    value = model(Object, var = var, x = timepoints)
    Object@master_tbl[,get_cols(Object, var, timepoints)] <- as.data.table(value)
  }
  Object
}


