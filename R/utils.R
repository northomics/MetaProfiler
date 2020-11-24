create_experimental_design <- function(filenames, dir, ...) {
  if(missing(filenames)) {
    filenames <- dir(dir, pattern = c("tsv|TSV|txt|TXT|csv|CSV|mzml|mzML|mzxml|mzXML"), full.names = T)
  }
  args <- list(...)
  d <- data.table::data.table(filenames = filenames, data.table::as.data.table(lapply(args, function(x) {
    stringi::stri_extract_last_regex(filenames, x)
  })))
  shush(d <- d[,lapply(.SD, function(x) {
    x_num <- as.numeric(x)
    if(all(is.na(x_num))) {
      return(x)
    }
    x_num
  })])
  d
}

plural <- function (singular, plural, n = NULL)
{
  return (ifelse(n>1, plural, singular))
}

.unequal_var_shrink <- function (stat, vars)
{
  n = length(stat)
  if (!n == length(vars))
    stop("stat and vars must be the same length.")
  vars[is.na(vars)] = 0
  s = stat[is.finite(stat) & is.finite(vars) & vars > 0]
  v = vars[is.finite(stat) & is.finite(vars) & vars > 0]
  mle_result = optimize(function(a_hat) {
    lik = -1 * (((s^2)/v) * (v/(v + a_hat))) + log(v/(v + a_hat))/2
    sum(c(1, lik))
  }, interval = c(0, 1e+06), maximum = TRUE)
  vars/(vars + mle_result$maximum)
}

.shrink <- function(data, method, vars, col, by, tu, trace) {
  by_t = c(tu, by)
  X <- data.table::copy(data)
  vars_trans = paste0(vars, "_trans")
  vars_trans_var = paste0(vars_trans, "_var")
  if(method == "ws") {
    X <- X[, c(vars_trans) := lapply(.SD, function(x) {
      test = .check_if_percantage_or_fraction(x)
      if(length(test) && test) {
        x = x/100
      }
      x = log(x/(1-x))
      odds = (1 - LFDR)/LFDR
      sw <- sum(odds, na.rm = T)
      xbar <- sum(odds * x, na.rm = T)/sw
      xvar = sum(odds * ((x - xbar)^2), na.rm = T)/(sw - 1)
      (x - xbar)
    }), by = c(tu), .SDcols = vars]
    X_var <- data.table::setnames(X[, lapply(.SD, function(x) {
      odds = (1 - LFDR)/LFDR
      sw <- sum(odds, na.rm = T)
      xbar <- sum(odds * x, na.rm = T)/sw
      sum(odds * ((x - xbar)^2), na.rm = T)/(sw - 1)
    }), by = by_t, .SDcols = vars_trans], c(by_t, vars_trans_var))
    X_mean <- X[, lapply(.SD, function(x) {
      odds = (1 - LFDR)/LFDR
      sw <- sum(odds, na.rm = T)
      sum(odds * x, na.rm = T)/sw
    }), by = by_t, .SDcols = c(col, vars_trans)]
    X <- X_var[X_mean, , on = c(by_t)]
    X = X[do.call(order, as.list(mget(by_t)))]
  } else {
    X <- X[, c(vars_trans) := lapply(.SD, function(x) {
      test = .check_if_percantage_or_fraction(x)
      if(length(test) && test) {
        x = x/100
      }
      x = log(x/(1-x))
      x_mean = mean(x, na.rm = T)
      (x - x_mean) # center the data
    }), by = tu, .SDcols = vars]
    X_var <- data.table::setnames(X[, lapply(.SD, var, na.rm = T), by = by_t, .SDcols = vars_trans],
                                 c(by_t, vars_trans_var))
    X_mean = X[, lapply(.SD, mean, na.rm = T), by = by_t, .SDcols = c(col, vars_trans)]
    X = X_var[X_mean, , on = c(by_t)]
    X = X[do.call(order, as.list(mget(by_t)))]
  }
  B_col = paste0(vars, " (B_hat)")
  B <- X[,c(by_t, vars_trans_var, vars_trans), with = F]
  B[,c(tu, B_col)] <- B[, lapply(vars, function(x) {
           x_trans = get(paste0(x, "_trans"))
           x_trans_var = get(paste0(x, "_trans_var"))
           .unequal_var_shrink(x_trans, x_trans_var)
         }), by = c(tu)]
  if(trace) {
    cat(crayon::blue("The average James-Stein shrinkage factor using unequal variance for ",
                     combine_words(vars), " is ", 
                     combine_words(signif(colMeans(B[,paste0(vars, " (B_hat)"), with = F]), 3), before = ""),
                     plural(".", ", respectively.", length(vars)), sep = ""))
  }
  X <- X[,c(by_t, col), with = F][B[,c(by_t, B_col),with=F], , on = c(by_t)]
  X <- X[, c(vars) := lapply(vars, function(x) {
      B = get(paste0(x, " (B_hat)"))
      # x_trans_var[is.na(x_trans_var)] = 0
      x_mean = get(x)
      ((1 - B) * x_mean + B * mean(x_mean, na.rm = T))
    })][,-..B_col]
  X
}


shrink <- function(Object, data = NULL, method, vars, by, timepoints = Object@timepoints, trace) {
  X = data.table::copy(Object@master_tbl)
  if(length(data)) {
    X = data.table::copy(data)
  }
  col = get_cols(Object, vars, x = timepoints)
  if(!is.list(timepoints)) {
    timepoints = list(timepoints)
  }
  if(length(timepoints) == 1 & length(vars) > 1) {
    timepoints = rep(timepoints, length(vars))
  }
  col_LFDR = get_cols(Object, rep("LFDR", length(timepoints)), x = timepoints)
  col_trans = paste0(col, "_trans")
  col_trans_var = paste0(col_trans, "_var")
  if(method == "ws") {
    X <- X[, c(col_trans) := mapply(function(x, c) {
      test = .check_if_percantage_or_fraction(x)
      if(length(test) && test) {
        x = x/100
      }
      x = log(x/(1-x))
      LFDR = get(c)
      odds = (1 - LFDR)/LFDR
      sw <- sum(odds, na.rm = T)
      xbar <- sum(odds * x, na.rm = T)/sw
      xvar = sum(odds * ((x - xbar)^2), na.rm = T)/(sw - 1)
      (x - xbar)/sqrt(xvar)
    }, .SD, col_LFDR, SIMPLIFY = F), .SDcols = col]
    X_var <- data.table::setnames(X[, mapply(function(x, c) {
      LFDR = get(c)
      odds = (1 - LFDR)/LFDR
      sw <- sum(odds, na.rm = T)
      xbar <- sum(odds * x, na.rm = T)/sw
      sum(odds * ((x - xbar)^2), na.rm = T)/(sw - 1)
    }, .SD, col_LFDR, SIMPLIFY = F), by = by, .SDcols = col_trans], c(by, col_trans_var))
    X_mean <- X[,mapply(function(x, c) {
      LFDR = get(c)
      odds = (1 - LFDR)/LFDR
      sw <- sum(odds, na.rm = T)
      sum(odds * x, na.rm = T)/sw
    }, .SD, col_LFDR, SIMPLIFY = F), by = by, .SDcols = c(col, col_trans)]
    X <- X_var[X_mean, , on = c(by)]
  } else {
    X <- X[, c(col_trans) := lapply(.SD, function(x) {
      test = .check_if_percantage_or_fraction(x)
      if(length(test) && test) {
        x = x/100
      }
      x = log(x/(1-x))
      x_mean = mean(x, na.rm = T)
      (x - x_mean) # center the data
    }), .SDcols = col]
    X_var = data.table::setnames(X[, lapply(.SD, var, na.rm = T), by = by, .SDcols = col_trans],
                                 c(by, col_trans_var))
    X_mean = X[, lapply(.SD, mean, na.rm = T), by = by, .SDcols = c(col, col_trans)]
    X = X_var[X_mean, , on = c(by)]
  }
  B_col = paste0(col, " (B_hat)")
  B <- X[,c(by, col_trans_var, col_trans), with = F]
  B[,B_col] <- B[, lapply(col, function(x) {
    x_trans = get(paste0(x, "_trans"))
    x_trans_var = get(paste0(x, "_trans_var"))
    .unequal_var_shrink(x_trans, x_trans_var)
  })]
  if(trace) {
    cat(crayon::blue("The average James-Stein shrinkage factor using unequal variance for ",
                     combine_words(col), " is ", 
                     combine_words(signif(colMeans(B[,paste0(col, " (B_hat)"), with = F]), 3), before = ""),
                     plural(".", ", respectively.", length(col)), "\n", sep = ""))
  }
  X <- X[,c(by, col), with = F][B[, c(by, B_col), with = F], , on = c(by)]
  X <- X[, c(col) := lapply(col, function(x) {
    B = get(paste0(x, " (B_hat)"))
    # x_trans_var[is.na(x_trans_var)] = 0
    x_mean = get(x)
    ((1 - B) * x_mean + B * mean(x_mean, na.rm = T))
  })][,-..B_col]
}

aggregate_by = function(Object, data = NULL, vars, by, timepoints = Object@timepoints,
                        method = c("weighted", "shrink", "ws", "mean", "median"), trace = T) {
  if(!length(data)) {
    data = Object@master_tbl
  }
  col = get_cols(Object, vars, x = timepoints)
  method = match.arg(method)
  dt <- data.table::copy(data)
  if(method == "weighted") {
    col_LFDR = get_cols(Object, rep("LFDR", length(timepoints)), x = timepoints)
    dt <- dt[, mapply(function(x, c) {
      LFDR = get(c)
      odds = (1 - LFDR)/LFDR
      sum(x * odds, na.rm = T)/sum(odds, na.rm = T)
    }, .SD, col_LFDR, SIMPLIFY = F), by = by, .SDcols = col]
  } else if (method == "shrink" | method == "ws") {
    dt <- shrink(Object, data, method, vars, by, timepoints, trace)
  } else if (method == "mean") {
    dt <- dt[, lapply(.SD, mean, na.rm = T), by = by, .SDcols = col]
  } else {
    dt <- dt[, lapply(.SD, median, na.rm = T), by = by, .SDcols = col]    
  }
  
  
}

.aggregate_by = function(data, by, tu, vars, col, method = c("weighted", "shrink", "ws", "mean", "median"), trace) {
  method = match.arg(method)
  dt <- data.table::copy(data)
  if(method == "weighted") {
    dt <- dt[, c(lapply(.SD, function(x) {
      LFDR[is.na(LFDR) & !is.na(x)] = 1
      odds = (1 - LFDR)/LFDR
      sum(x * odds)/sum(odds)
    }), list(LFDR = mean(LFDR, na.rm = T))), by = c(by, tu), .SDcols = col[-length(col)]]
  } else if (method == "shrink" | method == "ws") {
    dt <- .shrink(dt, method, vars, col, by, tu, trace)
  } else if (method == "mean") {
    dt <- dt[, lapply(.SD, mean, na.rm = T), by = c(by, tu), .SDcols = col]
  } else {
    dt <- dt[, lapply(.SD, median, na.rm = T), by = c(by, tu), .SDcols = col]    
  }
  dt
}



combine_words <- function (words, sep = ", ", and = " and ", before = "`", after = before, raw = T) 
{
  n = length(words)
  rs = xfun::raw_string
  if (n == 0) 
    return(words)
  words = paste0(before, words, after)
  if (n == 1) 
    return(rs(words))
  if (n == 2) 
    return(rs(paste(words, collapse = and)))
  if (grepl("^ ", and) && grepl(" $", sep)) 
    and = gsub("^ ", "", and)
  words[n] = paste0(and, words[n])
  if(raw)
   return(rs(paste(words, collapse = sep)))
  paste(words, collapse = sep)
}

rmsd_widths = function(Object, min_nb_groups = 4, min_nb_obs = 2){
  times = Object@timepoints
  windows = NULL
  scores = NULL
  for(i in 1:max(Object@timepoints)) {
    grp = NULL
    score = NULL
    for(t in Object@timepoints) {
      tmp = times[times >= t  & times <= (t + i)]
      if(!length(grp) || !any(sapply(grp, function(x) all(tmp %in% x)))) {
        grp = c(grp, list(tmp))
        score = c(score, diff(range(tmp)))
      }
    }
    score = sqrt(mean((score - mean(score))^2))
    if(!(score %in% scores) && length(grp) >= min_nb_groups && all(lengths(grp) >= min_nb_obs)) {
      windows = setNames(c(windows, list(grp)), c(names(windows), as.character(i)))
      scores = c(scores, score)
    }
  }
  list(scores = scores, windows = windows)
}

get_ranks <- function(x, main_ranks_only = T) {
  ranks = c("superkingdom", "kingdom", "subkingdom",
            "superphylum", "phylum", "subphylum",
            "superclass", "class", "subclass", "infraclass",
            "superorder", "order", "suborder", "infraorder", "parvorder",
            "superfamily", "family", "subfamily",
            "tribe", "subtribe",
            "genus", "subgenus",
            "species group", "species subgroup", "species", "subspecies",
            "varietas", "forma")
  if(main_ranks_only) {
    ranks = c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
    
  }
  colnames(x)[tolower(colnames(x)) %in% ranks]
}

shush <- function(expr, all = TRUE) {
  if (Sys.info()['sysname'] == "Windows") {
    file <- "NUL"
  } else {
    file <- "/dev/null"
  }
  
  if (all) {
    suppressWarnings(suppressMessages(suppressPackageStartupMessages(
      utils::capture.output(expr, file = file) 
    )))
  } else {
    utils::capture.output(expr, file = file)
  }
}

rescale <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)) {
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

save_plot <- function (filename, plot = last_plot(), device = NULL, path = NULL, 
                       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
                       dpi = 300, limitsize = TRUE, ...) {
  dpi <- ggplot2:::parse_dpi(dpi)
  dev <- ggplot2:::plot_dev(device, filename, dpi = dpi)
  dim <- ggplot2:::plot_dim(c(width, height), scale = scale, units = units, 
                            limitsize = limitsize)
  if (length(path)) {
    filename <- file.path(path, filename)
  }
  old_dev <- grDevices::dev.cur()
  dev(filename = filename, width = dim[1], height = dim[2])
  on.exit(utils::capture.output({
    grDevices::dev.off()
    if (old_dev > 1) grDevices::dev.set(old_dev)
  }))
  if(inherits(plot, c("Heatmap", "AnnotationFunction", "Legends", "HeatmapAnnotation", "SingleAnnotation", "HeatmapList"))) {
    ComplexHeatmap::draw(plot, ...)
  } else if (inherits(plot, "igraph")) {
    plot(plot, ...)
  } else if (inherits(plot, "gtable")) {
    gridExtra::grid.arrange(plot, ...)
  } else {
    grid::grid.draw(plot)
  }
  invisible()
}

.guess_accession_column <- function(d, trace, pattern) {
  cols <- toupper(colnames(d))
  suppressWarnings(idx <- which(grepl("ACC|PRO|QUERY|NAME|GENE|LEAD|RAZOR|GROUP", cols, perl = T)))
  if(length(idx) == 0) {
    cat(crayon::red(paste0("Could not guess the name of the accession/protein column for `", deparse(substitute(d)), "`. Please provide the name of this column.\n")))
    return(NULL)
  }
  lookahead <- NULL
  for(i in idx) {
    c = d[,get(colnames(d)[i])]
    if(is.character(c)) {
      lookahead <- c(lookahead, sum(stringi::stri_count_regex(d[[colnames(d)[i]]], pattern)/nchar(d[[colnames(d)[i]]]), na.rm = T))
    } else {
      lookahead <- c(lookahead, -Inf)
    }
  }
  if(max(lookahead, na.rm = TRUE) > 0) {
    idx <- idx[which.max(lookahead)]
    if(trace) {
      cat(crayon::silver(paste0("Using `", colnames(d)[idx], "` as the accession/protein column.\n")))
    }
    return(colnames(d)[idx])
  } else {
    cat(crayon::red(paste0("Could not guess the name of the accession/protein column for `", deparse(substitute(d)), "`. Please provide the name of this column.\n")))
    return(NULL)
  }
}

.convert_modified_peptide_column <- function(d, col, trace) {
  if(trace) {
    cat(crayon::blue(paste0("Looks like `", col, "` is a peptide column, but appears to contain post translational modifications.\n")))
  }
  new_column = paste(col, "(no PTMs)")
  if(trace) {
    cat(crayon::blue(paste0("Creating column `", new_column,"` with sequences from column `", col, "`.\n")))
  }
  d[,new_column] = d[[col]]
  if(trace) {
    cat(crayon::blue(paste0("Removing PTMs...")))
  }
  convert <- unique(d[,new_column,with=F])
  convert[,"seq"] <- .remove_modifications(convert[[new_column]])
  
  if(all(!grepl("[^ILKMFTWVRHANDCEQGPSY]", convert$seq, perl = T)))
  {
    d[,new_column] <- convert[d, seq, on = new_column]
  } else {
    cat(crayon::black(paste0("Failed to remove PTMs, returning NULL.\n")))
    return(NULL)
  }
  if(trace) {
    cat(crayon::blue("done.\n"))
  }
  list(d, new_column)
}

.guess_unmodified_peptide_column <- function(d, trace) {
  cols <- toupper(colnames(d))
  suppressWarnings(idx <- which(grepl("PEP|SEQ", cols)))
  if(!all(is.finite(idx))) {
    warning(paste0("Could not guess the name of the unmodified peptide column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
  lookahead <- NULL
  for(i in idx) {
    c = d[,get(colnames(d)[i])]
    if(is.character(c)) {
      lookahead <- c(lookahead, sum(!grepl("[^ILKMFTWVRHANDCEQGPSY]", c, perl = T)))
    } else {
      lookahead <- c(lookahead, 0)
    }
  }
  idx2 <- idx[lookahead == nrow(d)]
  if(length(idx2) > 0) {
    idx <- idx[which.max(lookahead == nrow(d))]
    if(trace) {
      cat(crayon::silver(paste0("Using `", colnames(d)[idx], "` as the unmodified peptide column.\n")))
    }
    return(colnames(d)[idx])
  } else {
    idx2 <- idx[lookahead > (nrow(d) * 0.5)]
    if(length(idx2) > 0) {
      idx <- idx[which.max(lookahead > (nrow(d) * 0.5))]
      res <- .convert_modified_peptide_column(d, colnames(d)[idx], trace)
      if (trace) {
        cat(crayon::silver(paste0("Using `", res[[2]], "` as the unmodified peptide column.\n")))
      }
      return(res)
    }
    warning(paste0("Could not guess the name of the unmodified peptide column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
}

.guess_modified_peptide_column <- function(d, trace) {
  cols <- toupper(colnames(d))
  suppressWarnings(idx <- which(grepl("PEP|SEQ", cols)))
  if(!all(is.finite(idx))) {
    warning(paste0("Could not guess the name of the modified peptide column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
  
  lookahead <- NULL
  for(i in idx) {
    c = d[,get(colnames(d)[i])]
    if(is.character(c)) {
      c2 <- .remove_modifications(c)
      lookahead <- c(lookahead, list(c(sum(!grepl("[^ILKMFTWVRHANDCEQGPSY]", c, perl = T)), sum(grepl("[\\(\\)\\.\\{\\}\\<\\>\\[\\]\\_]", c, perl = T)))))
    } else {
      lookahead <- c(lookahead,list(c(-Inf, -Inf)))
    }
  }
  idx <- idx[sapply(lookahead, function(x) (x[1] > 0) && (x[2] > 0))]
  if(length(idx) > 0) {
    if(trace) {
      cat(crayon::silver(paste0("Using `", colnames(d)[idx][1], "` as the modified peptide column.\n")))
    }
    colnames(d)[idx][1]
  } else {
    warning(paste0("Could not guess the name of the modified peptide column for `", deparse(substitute(d)), "`. Please provide the name of this column."))
    return(NULL)
  }
}

.get_columns <- function(x, cols) {
  x = gsub('([\\[!@#$%^&*(),.?":{}|<>\\]])', "\\\\\\1", x, perl = T)
  cols[grepl(paste0("^", x,"\\s?\\S+"), cols)]
}

guess_measurements <- function(x, cols) {
  if(length(x) == 0 || is.na(x)) return(NULL)
  cols <- .get_columns(x, cols)
  if(length(unique(cols)) == 0)
  {
    warning(paste0("No valid columns that contains the variable `", x,"`. Please add columns where the first word is the variable's name, followed by a unique identifier (e.g. ",
                   paste0(x, " ", 1:3, ", ", collapse = ""), "..., or ", paste0(x, " ", c("light", "heavy"), ", ", collapse = ""),"). Returning NULL."))
    return(NULL)
  }
  cols
}

.check_columns <- function(x, cols) {
  if(!length(x) || is.na(x)) return(character())
  cols <- .get_columns(x, cols)
  # if(length(unique(cols)) == 0)
  # {
    # stop(paste0("No valid columns that contains the variable `", x,"`. Please add columns where the first word is the variable's name, followed by a unique identifier (e.g. ",
                # paste0(x, " ", 1:3, ", ", collapse = ""), "..., or ", paste0(x, " ", c("light", "heavy"), ", ", collapse = ""),")."))
    
  # }
  cols
}

# Removes PTMs by deleting characters within brackets and by deleting every character that is not a one-letter code for an amino acid...
.remove_modifications <- function(x) {
  stringi::stri_replace_all_regex(stringi::stri_replace_all_regex(x, "\\([\\S\\s]+?\\)|\\[[\\S\\s]+?\\]|\\{[\\S\\s]+?\\}|\\<[\\S\\s]+?\\>|_[\\S\\s]+?_", ""),
                                  "[^ILKMFTWVRHANDCEQGPSY]", "")
}

.check_if_percantage_or_fraction <- function(x) {
  shush(v <- max(x, na.rm = T))
  if(v > 1 && v <= 100) {
    return(T)
  } else if(v > 0 && v <= 1) {
    return(F)
  } else {
    return(NULL)
  }
}

.index_column <- function(data, idx, default, trace) {
  col = colnames(data)[idx]
  if(!length(col)) {
    if(trace) {
      cat(crayon::blue(paste0("There is no column name at ", deparse(substitute(data)), "[", idx, "]. Renaming this column as ", default, ".")))
    }
    return(default)
  }
  return(col)
}



.guess_columns <- function(
  data,
  data_peptide_column_PTMs,
  data_peptide_column_no_PTMs,
  data_accession_column,
  accession_pattern,
  annotate_by_peptide_with_PTMs,
  
  pep2pro,
  pep2pro_peptide_column,
  pep2pro_accession_column,
  pep2pro_accession_pattern,
  
  pep2func,
  pep2func_peptide_column,
  pep2func_function_columns,
  
  pep2taxon,
  pep2taxon_peptide_column,
  rank_columns,
  main_ranks_only,
  
  pro2func,
  pro2func_accession_column,
  pro2func_accession_pattern,
  pro2func_function_columns,
  
  trace
) {
  
  
  if(!missing(data) && length(data)){
    if(any(c(data_peptide_column_PTMs,
             data_peptide_column_no_PTMs,
             data_accession_column) == "guess") && trace) {
      cat(crayon::red("***** Protein-SIP Table *****\n"))
    }
    if(length(data_peptide_column_PTMs) && data_peptide_column_PTMs == "guess") {
      data_peptide_column_PTMs <- .guess_modified_peptide_column(data, trace)
    } else if (is.numeric(data_peptide_column_PTMs)) {
      tmp = .index_column(data, data_peptide_column_PTMs, "Peptides (PTMs)", trace)
      colnames(data)[data_peptide_column_PTMs] = tmp
      data_peptide_column_PTMs = tmp
    }
    if (length(data_peptide_column_no_PTMs) && data_peptide_column_no_PTMs == "guess") {
      data_peptide_column_no_PTMs = .guess_unmodified_peptide_column(data, trace)
      if(is.list(data_peptide_column_no_PTMs)) {
        data = data_peptide_column_no_PTMs[[1]]
        data_peptide_column_no_PTMs = data_peptide_column_no_PTMs[[2]]
      }
    } else if (is.numeric(data_peptide_column_no_PTMs)) {
      tmp = .index_column(data, data_peptide_column_no_PTMs, "Peptides (no PTMs)", trace)
      colnames(data)[data_peptide_column_no_PTMs] = tmp
      data_peptide_column_no_PTMs = tmp
    }
    # if(!length(data_peptide_column_PTMs) && !length(data_peptide_column_no_PTMs)) {
      # stop("Both variable `data_peptide_column_PTMs` and `data_peptide_column_no_PTMs` cannot be NULL")
    # }
    if (!length(data_peptide_column_PTMs) && length(data_peptide_column_no_PTMs)) {
      data_peptide_column_PTMs <- paste(data_peptide_column_no_PTMs, "(PTMs)")
      if(trace) {
        cat(crayon::blue(paste0("Copying and pasting sequences in column `", data_peptide_column_no_PTMs, 
                                "` to column `", data_peptide_column_PTMs, "`.\n")))
      }
      data[[data_peptide_column_PTMs]] = data[[data_peptide_column_no_PTMs]]
    }
    if (length(data_peptide_column_no_PTMs) && !length(data_peptide_column_no_PTMs)) {
      data_peptide_column_no_PTMs <- paste(data_peptide_column_PTMs, "(no PTMs)")
      if(trace) {
        cat(crayon::blue(paste0("Copying and pasting sequences in column `", data_peptide_column_PTMs, "` to column `", data_peptide_column_no_PTMs, "`.\n")))
      }
      data[,data_peptide_column_no_PTMs] = data[[data_peptide_column_PTMs]]
      if(trace) {
        cat(crayon::blue(paste0("Removing PTMs in column `", data_peptide_column_no_PTMs, "`...")))
      }
      convert <- unique(data[,data_peptide_column_no_PTMs,with=F])
      convert[,"seq"] <- .remove_modifications(convert[[data_peptide_column_no_PTMs]])
      data[,data_peptide_column_no_PTMs] <- convert[data, seq, on = data_peptide_column_no_PTMs]
      if(trace) {
        cat(crayon::blue("done.\n"))
      }
    }
    
    annotate_with = ifelse(
      !annotate_by_peptide_with_PTMs,
      data_peptide_column_no_PTMs,
      data_peptide_column_PTMs
    )
    
    rows_with_PTMs = grepl("[^ILKMFTWVRHANDCEQGPSY]", data[[data_peptide_column_no_PTMs]])
    if(any(rows_with_PTMs)) {
      res = .convert_modified_peptide_column(data, data_peptide_column_no_PTMs, trace)
      if(is.list(res)) {
        data = res[[1]]
        data_peptide_column_no_PTMs = res[[2]]
      }
    }
    if(length(data_accession_column) && data_accession_column == "guess") {
      data_accession_column <- .guess_accession_column(data, trace, accession_pattern)
    } else if (is.numeric(data_accession_column)) {
      tmp = .index_column(data, data_accession_column, "Proteins", trace)
      colnames(data)[data_accession_column] = tmp
      data_accession_column = tmp
    }
    # 
    # if(!length(data_accession_column)) {
    #   data_accession_column <- "Proteins"
    # }
  } else {
    data_peptide_column_PTMs = NULL
    data_peptide_column_no_PTMs = NULL
    data_accession_column = NULL
    accession_pattern = NULL
    annotate_by_peptide_with_PTMs = NULL
  }
  
  
  
  if(!missing(pep2pro) && length(pep2pro)) {
    
    if(!inherits(pep2pro, "data.frame")) {
      stopifnot(file.exists(pep2pro))
      pep2pro <- data.table::fread(pep2pro)
    }
    
    if(any(c(pep2pro_peptide_column,
             pep2pro_accession_column) == "guess") && trace) {
      cat(crayon::red("***** Peptide-to-Protein Table *****\n"))
    }
    if(length(pep2pro_accession_column) && pep2pro_accession_column == "guess") {
      pep2pro_accession_column <- .guess_accession_column(pep2pro, trace, pep2pro_accession_pattern)
    } else if (is.numeric(pep2pro_accession_column)) {
      tmp = .index_column(pep2pro, pep2pro_accession_column, "Proteins", trace)
      colnames(pep2pro)[pep2pro_accession_column] = tmp
      pep2pro_accession_column = tmp
    }
    
    if(annotate_with == data_peptide_column_no_PTMs){
      if(pep2pro_peptide_column == "guess") {
        pep2pro_peptide_column <- .guess_unmodified_peptide_column(pep2pro, trace)
        if(is.list(pep2pro_peptide_column)) {
          pep2pro = pep2pro_peptide_column[[1]]
          pep2pro_peptide_column = pep2pro_peptide_column[[2]]
        }
      } else if (is.numeric(pep2pro_peptide_column)) {
        tmp = .index_column(pep2pro, pep2pro_peptide_column, "Peptides (no PTMs)", trace)
        colnames(pep2pro)[pep2pro_peptide_column] = tmp
        pep2pro_peptide_column = tmp
      }
      rows_with_PTMs = grepl("[^ILKMFTWVRHANDCEQGPSY]", pep2pro[[pep2pro_peptide_column]])
      if(any(rows_with_PTMs)) {
        res = .convert_modified_peptide_column(pep2pro, pep2pro_peptide_column, trace)
        if(is.list(res)) {
          pep2pro = res[[1]]
          pep2pro_peptide_column = res[[2]]
        }
      }
    } else if(annotate_with == data_peptide_column_PTMs && pep2pro_peptide_column == "guess") {
      pep2pro_peptide_column <- .guess_modified_peptide_column(pep2pro, trace)
    } else if (is.numeric(pep2pro_peptide_column)) {
      tmp = .index_column(pep2pro, pep2pro_peptide_column, "Peptides (PTMs)", trace)
      colnames(pep2pro)[pep2pro_peptide_column] = tmp
      pep2pro_peptide_column = tmp
    }
    rows_with_aa = grepl("[ILKMFTWVRHANDCEQGPSY]", pep2pro[[pep2pro_peptide_column]])
    if(!all(rows_with_aa)) {
      warning(paste0("Row(s) [", ifelse(sum(!rows_with_aa) > 5, 
                                     paste0(paste0(which(!rows_with_aa)[1:5], collapse = ","), "..."),
                                     paste0(which(!rows_with_aa), collapse = ",")),
                  "] in column `", pep2pro_peptide_column, "` for `pep2pro` seem to contain characters that are not single letter amino acid codes"))
    }
  } else {
    pep2pro_peptide_column = NULL
    pep2pro_accession_column = NULL
  }
  
  
  if(!missing(pep2taxon) && length(pep2taxon)) {
    
    if(!inherits(pep2taxon, "data.frame")) {
      print(pep2taxon)
      stopifnot(file.exists(pep2taxon))
      pep2taxon <- data.table::fread(pep2taxon)
    }
    
    if(any(c(pep2taxon_peptide_column,
             rank_columns) == "guess") && trace) {
      cat(crayon::red("***** Peptide-to-Taxon Table *****\n"))
    }
    if(annotate_with == data_peptide_column_no_PTMs){
      if(pep2taxon_peptide_column == "guess") {
        pep2taxon_peptide_column <- .guess_unmodified_peptide_column(pep2taxon, trace)
        if(is.list(pep2taxon_peptide_column)) {
          pep2taxon = pep2taxon_peptide_column[[1]]
          pep2taxon_peptide_column = pep2taxon_peptide_column[[2]]
        }
      } else if (is.numeric(pep2taxon_peptide_column)) {
        tmp = .index_column(pep2taxon, pep2taxon_peptide_column, "Peptides (no PTMs)", trace)
        colnames(pep2taxon)[pep2taxon_peptide_column] = tmp
        pep2taxon_peptide_column = tmp
      }
      rows_with_PTMs = grepl("[^ILKMFTWVRHANDCEQGPSY]", pep2taxon[[pep2taxon_peptide_column]])
      if(any(rows_with_PTMs)) {
        res = .convert_modified_peptide_column(pep2taxon, pep2taxon_peptide_column, trace)
        if(is.list(res)) {
          pep2taxon = res[[1]]
          pep2taxon_peptide_column = res[[2]]
        }
      }
      rows_with_aa = grepl("[ILKMFTWVRHANDCEQGPSY]", pep2taxon[[pep2taxon_peptide_column]])
      if(!all(rows_with_aa)) {
        warning(paste0("Row(s) [", ifelse(
          sum(!rows_with_aa) > 5, 
          paste0(paste0(which(!rows_with_aa)[1:5], collapse = ","), "..."),
          paste0(which(!rows_with_aa), collapse = ",")),
          "] in column `", pep2taxon_peptide_column, 
          "` for `pep2taxon` seem to contain characters that are not single letter amino acid codes"))
      }
    } else if(annotate_with == data_peptide_column_PTMs && pep2taxon_peptide_column == "guess") {
      pep2taxon_peptide_column <- .guess_modified_peptide_column(pep2taxon, trace)
    } else if (is.numeric(pep2taxon_peptide_column)) {
      tmp = .index_column(pep2taxon, pep2taxon_peptide_column, "Peptides (PTMs)", trace)
      colnames(pep2taxon)[pep2taxon_peptide_column] = tmp
      pep2taxon_peptide_column = tmp
    }
    if((length(rank_columns) == 1) && (rank_columns == "guess")) {
      rank_columns = get_ranks(pep2taxon, main_ranks_only)
    }
    pep2taxon_columns = colnames(pep2taxon)
    lca_column = pep2taxon_columns[grepl("lca", tolower(pep2taxon_columns))]
    lca_rank_column = pep2taxon_columns[grepl("rank", tolower(pep2taxon_columns))]
    pep2taxon[pep2taxon == ""] <- NA
    if(length(lca_rank_column) == 0) {
      pep2taxon <- pep2taxon[, rank := colnames(.SD)[max(which(!is.na(unlist(.SD))))], by = 1:nrow(pep2taxon), .SDcols = rank_columns]
      lca_rank_column = "rank"
    }
    if(length(lca_column) == 0) {
      pep2taxon <- pep2taxon[, lca := .SD[, max(which(!is.na(unlist(.SD)))), with = F], by = 1:nrow(pep2taxon), .SDcols = rank_columns]
      lca_column = "lca"
    }
    pep2taxon_columns <- c(lca_column, lca_rank_column, rank_columns)
  } else {
    pep2taxon_columns = NULL
    rank_columns = NULL
  }

  
  if(!missing(pep2func) && length(pep2func)) {
    
    if(!inherits(pep2func, "data.frame")) {
      stopifnot(file.exists(pep2func))
      pep2func <- data.table::fread(pep2func)
    }
    
    if(any(c(pep2func_peptide_column, pep2func_function_columns) == "guess") && trace) {
      cat(crayon::red("***** Peptide-to-Function Table *****\n"))
    }
    if((length(pep2func_function_columns) == 1) && 
       (pep2func_function_columns == "guess")) {
      pep2func_function_columns <- colnames(pep2func)[grepl("^COG|^NOG|^KEGG|^GO|^BRITE|^REACTOME|^PRO[\\s\\S]+?NAME|^INTERPRO|^EC", toupper(colnames(pep2func)), perl = T)]
      if(trace) {
        cat(crayon::silver(paste0("Using ", combine_words(pep2func_function_columns)," as the functional annotation ",
                                  plural("column", "columns", length(pep2func_function_columns)),"\n")))
      }
    }
    if(annotate_with == data_peptide_column_no_PTMs){
      if(pep2func_peptide_column == "guess") {
        pep2func_peptide_column <- .guess_unmodified_peptide_column(pep2func, trace)
        if(is.list(pep2func_peptide_column)) {
          pep2func = pep2func_peptide_column[[1]]
          pep2func_peptide_column = pep2func_peptide_column[[2]]
        }
      } else if (is.numeric(pep2func_peptide_column)) {
        tmp = .index_column(pep2func, pep2func_peptide_column, "Peptides (no PTMs)", trace)
        colnames(pep2func)[pep2func_peptide_column] = tmp
        pep2func_peptide_column = tmp
      }
      rows_with_PTMs = grepl("[^ILKMFTWVRHANDCEQGPSY]", pep2func[[pep2func_peptide_column]])
      if(any(rows_with_PTMs)) {
        res = .convert_modified_peptide_column(pep2func, pep2func_peptide_column, trace)
        if(is.list(res)) {
          pep2func = res[[1]]
          pep2func_peptide_column = res[[2]]
        }
      }
      rows_with_aa = grepl("[ILKMFTWVRHANDCEQGPSY]", pep2func[[pep2func_peptide_column]])
      if(!all(rows_with_aa)) {
        warning(paste0("Row(s) [", ifelse(
          sum(!rows_with_aa) > 5, 
          paste0(paste0(which(!rows_with_aa)[1:5], collapse = ","), "..."),
          paste0(which(!rows_with_aa), collapse = ",")),
          "] in column `", pep2func_peptide_column, 
          "` for `pep2func` seem to contain characters that are not single letter amino acid codes"))
      }
    } else if(annotate_with == data_peptide_column_PTMs && pep2func_peptide_column == "guess") {
      pep2func_peptide_column <- .guess_modified_peptide_column(pep2func, trace)
    } else if (is.numeric(pep2func_peptide_column)) {
      tmp = .index_column(pep2func, pep2func_peptide_column, "Peptides (PTMs)", trace)
      colnames(pep2func)[pep2func_peptide_column] = tmp
      pep2func_peptide_column = tmp
    }
  } else {
    pep2func_peptide_column = NULL
    pep2func_function_columns = NULL
  }
  
  
  if(!missing(pro2func) && length(pro2func)) {
    
    if(!inherits(pro2func, "data.frame")) {
      stopifnot(file.exists(pro2func))
      pro2func <- data.table::fread(pro2func)
    }
    
    if(any(c(pro2func_accession_column, pro2func_function_columns) == "guess") && trace) {
      cat(crayon::red("***** Protein-to-Function Table *****\n"))
    }
    if((length(pro2func_function_columns) == 1) && 
       (pro2func_function_columns == "guess")) {
      pro2func_function_columns <- colnames(pro2func)[grepl("^COG|^NOG|^KEGG|^GO|^BRITE|^REACTOME|^PRO[\\s\\S]+?NAME", toupper(colnames(pro2func)), perl = T)]
      if(trace) {
        cat(crayon::silver(paste0("Using ", combine_words(pro2func_function_columns)," as the functional annotation ",
                                  plural("column", "columns", length(pro2func_function_columns)),"\n")))
      }
    }
    if(length(pro2func_accession_column) && (pro2func_accession_column == "guess")) {
      pro2func_accession_column <- .guess_accession_column(pro2func, trace, pro2func_accession_pattern)
    } else if (is.numeric(pro2func_accession_column)) {
      tmp = .index_column(pro2func, pro2func_accession_column, "Proteins", trace)
      colnames(pro2func)[pro2func_accession_column] = tmp
      pro2func_accession_column = tmp
    }
  } else {
    pro2func_accession_column = NULL
    pro2func_function_columns = NULL
  }
  
  list(data = data,
       data_peptide_column_PTMs = data_peptide_column_PTMs,
       data_peptide_column_no_PTMs = data_peptide_column_no_PTMs,
       data_accession_column = data_accession_column,
       annotate_with = annotate_with,
       pep2pro = pep2pro,
       pep2pro_peptide_column = pep2pro_peptide_column,
       pep2pro_accession_column = pep2pro_accession_column,
       pep2func = pep2func,
       pep2func_peptide_column = pep2func_peptide_column,
       pep2func_function_columns = pep2func_function_columns,
       pep2taxon = pep2taxon,
       pep2taxon_peptide_column = pep2taxon_peptide_column,
       rank_columns = rank_columns,
       pep2taxon_columns = pep2taxon_columns,
       pro2func = pro2func,
       pro2func_accession_column = pro2func_accession_column,
       pro2func_function_columns = pro2func_function_columns)
  
}

store <- function(mz,RT,sequence,charge,aa_before,aa_after,start,end,ids,scan,peptide_score,protein_score = 0,
                  pro2id, output, db, Experiment, document_id, used_target_decoy = T, charges = "1+, 4+",
                  fixed_modifications = c("Carbamidomethyl (C)"),
                  variable_modifications = c("Acetyl (N-term)", "Oxidation (M)"),
                  taxnomomy = "", enzyme = "trypsin",
                  mass_type = "monoisotopic", missed_cleavages = 2, precursor_peak_tolerance = 0,
                  precursor_peak_tolerance_ppm = F, peak_mass_tolerance = 0, peak_mass_tolerance_ppm = F,
                  peptide_score_type = "", peptide_higher_score_better = F,
                  protein_score_type = "", protein_higher_score_better= F, protein_charges = "")
{
  require(dplyr)
  os <- paste(
    # IdXML header
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",
    "<?xml-stylesheet type=\"text/xsl\" href=\"https://www.openms.de/xml-stylesheet/IdXML.xsl\" ?>\n",
    "<IdXML version=\"1.5\"",
    " id=\"", document_id, "\"",
    " xsi:noNamespaceSchemaLocation=\"https://www.openms.de/xml-schema/IdXML_1_5.xsd\" ",
    "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n",
    
    # IdXML search parameters.
    "\t<SearchParameters charges=\"", charges, "\"",
    " id=\"SP_0\"",
    " db=\"", paste0(db, collapse = " "), "\"",
    " db_version=\"0\"",
    " taxnomomy=\"", taxnomomy, "\"",
    " mass_type=\"", mass_type, "\"",
    " enzyme=\"", enzyme, "\"",
    " missed_cleavages=\"", missed_cleavages, "\"",
    " precursor_peak_tolerance=\"", precursor_peak_tolerance, "\"",
    " precursor_peak_tolerance_ppm=\"", ifelse(precursor_peak_tolerance_ppm, "true", "false"), "\"",
    " peak_mass_tolerance=\"", peak_mass_tolerance, "\"", "",
    " peak_mass_tolerance_ppm=\"", ifelse(peak_mass_tolerance_ppm, "true", "false"), "\" >\n",
    paste0("\t\t<FixedModification name=\"", fixed_modifications, "\" />", collapse = "\n"),
    "\n",
    paste0("\t\t<VariableModification name=\"", variable_modifications ,"\" />", collapse = "\n"),
    "\n",
    "\t\t\t<UserParam type=\"string\" name=\"TargetDecoyApproach\" value=\"",
    ifelse(used_target_decoy, "true", "false"), "\"/>\n",
    "\t</SearchParameters>\n",
    
    "\t<IdentificationRun",
    " date=\"", gsub(" ", "T", as.character(Sys.time())), "\"",
    " search_engine=\"MetaLab\"",
    " search_engine_version=\"2.0.0\"",
    " search_parameters_ref=\"SP_0\"",
    " >\n",
    
    "\t\t<ProteinIdentification",
    " score_type=\"", protein_score_type, "\"",
    " charges=\"0\"",
    " higher_score_better=\"", ifelse(protein_higher_score_better, "true", "false"), "\"",
    " significance_threshold=\"0\" >\n",
    paste0("\t\t\t<ProteinHit ",
           "id=\"PH_" ,  pro2id$id,
           "\" accession=\"", pro2id$protein,
           "\" score=\"", protein_score, "\" ",
           "sequence=\"", pro2id$sequence, "\" >\n",
           "\t\t\t\t<UserParam type=\"string\" name=\"target_decoy\" value=\"",
           case_when(grepl("DECOY_|REV_", pro2id$protein) ~ "decoy", T ~ "target"), "\"/>\n",
           "\t\t\t\t<UserParam type=\"string\" name=\"Description\" value=\"", pro2id$description, "\"/>\n",
           "\t\t\t</ProteinHit>",
           collapse = "\n"),
    "\n\t\t</ProteinIdentification>\n",
    
    
    paste0("\t\t<PeptideIdentification ",
           "score_type=\"", peptide_score_type, "\" ",
           "higher_score_better=\"", ifelse(peptide_higher_score_better, "true", "false"), "\" ",
           "significance_threshold=\"0\" ",
           "MZ=\"" ,  mz ,  "\" ",
           "RT=\"" ,  RT,  "\" >\n",
           "\t\t\t<PeptideHit",
           " score=\"", peptide_score, "\"",
           " sequence=\"" ,  sequence,  "\"",
           " charge=\"" , charge,  "\"",
           " aa_before=\"", aa_before, "\"",
           " aa_after=\"", aa_after, "\"",
           " start=\"", start, "\"",
           " end=\"", end, "\"",
           " protein_refs=\"", ids, "\"",
           " spectrum_reference=\"scan=",scan,"\"",
           " >\n",
           "\t\t\t</PeptideHit>\n",
           "\t\t</PeptideIdentification>",
           collapse = "\n"),
    
    "\n\t</IdentificationRun>\n",
    "</IdXML>\n", sep = ""
  )
  if(!dir.exists(output)) {
    dir.create(output)
  }
  writeLines(os, con = file.path(output, paste(Experiment, ".idXML", sep = "")))
}

MQ2idXML <- function(MQ.output, db, output)
{
  db = "D:/sample_specific_db_decoy.fasta"
  MQ.output = "C:/Users/Patrick Smyth/OneDrive - University of Ottawa/supplemental_tables/pipeline_output/metalab_output/combined_12_6/txt"
  peptides = data.table::fread(file.path(MQ.output, "peptides.txt"))
  peptides = peptides[peptides$`Potential contaminant` != "+" & peptides$Proteins != "",]
  
  msms = data.table::fread(file.path(MQ.output, "msms.txt"))
  msms$`Modified sequence` <- gsub("\\(ac\\)", "(Acetyl)", msms$`Modified sequence`)
  msms$`Modified sequence` <- gsub("\\(ox\\)", "(Oxidation)", msms$`Modified sequence`)
  msms$`Modified sequence` <- gsub("_\\(", ".(", msms$`Modified sequence`)
  msms$`Modified sequence` <- gsub("_", "", msms$`Modified sequence`)
  msms$`Modified sequence` <- gsub("C", "C(Carbamidomethyl)", msms$`Modified sequence`)
  msms = msms[Sequence %in%  peptides$Sequence]
  annotation_table <- unique(msms[, .(Sequence, Proteins)])
  annotation_table$Proteins <- stringi::stri_split_fixed(annotation_table$Proteins, ";")
  annotation_table <- annotation_table[, .(Proteins = unlist(Proteins)), by = Sequence]
  
  proteins <- unique(annotation_table$Proteins)
  sourceCpp("~MetaProfiler/fasta_reader.cpp")
  fasta <- data.table::as.data.table(get_accessions(db, proteins))
  fasta$sequence <- gsub("\n", "", fasta$sequence)
  annotation_table$`Protein Sequence` <- fasta[annotation_table$Proteins,sequence, on = "accession"]
  annotation_table$`Protein Descriptions` <- fasta[annotation_table$Proteins,description, on = "accession"]
  annotation_table <- na.omit(annotation_table, "Protein Sequence")
  pro_seq = annotation_table$`Protein Sequence`
  pep_seq = annotation_table$Sequence
  pos = stringi::stri_locate_all_fixed(
    pro_seq, pep_seq
  )
  chars = strsplit(annotation_table$`Protein Sequence`, "")
  add = mapply(function(p, s) {
    s = c("[", s, "]")
    start = p[,1]
    end = p[,2]
    b = start
    a = end + 2
    aa_before = paste0(s[b], collapse = " ")
    aa_after = paste0(s[a], collapse = " ")
    start = paste0(start - 1, collapse = " ")
    end = paste0(end - 1, collapse = " ")
    list(aa_before = aa_before, aa_after = aa_after, start = start, end = end)
  }, pos, chars, SIMPLIFY = F)
  add = data.table::rbindlist(add)
  test = rowSums(add != "NA") != 0
  annotation_table[, colnames(add)] = add
  annotation_table = annotation_table[test]
  annotation_table = unique(annotation_table)
  col.names <- colnames(annotation_table)
  annotation_table <- annotation_table[, c(lapply(.SD[, col.names[2:4], with = F], list), lapply(.SD[, col.names[5:8], with = F], paste0, collapse = " ")), by = Sequence]
  msms <- annotation_table[msms, , on = "Sequence"]
  msms[is.na(msms)] <- 0
  output = "D:/knime_out/idxmls"
  for (e in unique(msms$`Raw file`))
  {
    x = msms[`Raw file` == e,]
    pro2id <- data.table::data.table(protein = unlist(x$Proteins), sequence = unlist(x$`Protein Sequence`),
                         description = unlist(x$`Protein Descriptions`))
    pro2id <- unique(pro2id)
    pro2id$id <- 1:nrow(pro2id)
    map2pep <- rep(1:nrow(x), lengths(x$Proteins))
    ids <- pro2id[unlist(x$Proteins), id, on = "protein"]
    ids <- split(ids, map2pep)
    ids <- sapply(ids, function(x) paste0("PH_", x, collapse = " "))
    x$Protein_ids <- ids
    store(
      mz = x$`m/z`,
      RT = x$`Retention time` * 60,
      sequence = x$`Modified sequence`,
      charge = x$Charge,
      aa_before = x$aa_before,
      aa_after = x$aa_after,
      start = x$start,
      end = x$end,
      document_id = e,
      Experiment = e,
      ids = ids,
      scan = x$`Scan number`,
      peptide_score = x$PEP,
      peptide_score_type = "PEP",
      peptide_higher_score_better = T,
      pro2id = pro2id,
      output = output,
      db = db 
    )
    
  }
}

file = "D:/rhizosphere/soluble.csv"
db = "D:/db/db1.fasta"
output = "D:/rhizosphere/idxml"
Scaffold2idXML <- function(file, db, output, file_id = NULL, pep_time_0_only = T,
                           time_0_contains = "0h")
{
  csv = lapply(file, data.table::fread)
  if(length(file_id))
    csv = mapply(function(x, y) {
      x$`Biological sample category` = paste0(y, "_", x$`Biological sample category`)
      x$`Biological sample name` = paste0(y, "_", x$`Biological sample name`)
      x[-nrow(x)]
    }, csv, file_id, SIMPLIFY = F)
  csv = data.table::rbindlist(csv, fill = T)
  csv[csv == ""] <- NA
  csv$Sequence = Sequence = gsub("m", "M", csv$`Peptide sequence`)
  csv$`Peptide sequence` <- gsub("m", "M(Oxidation)", csv$`Peptide sequence`)
  if(pep_time_0_only) {
    seq_0 = csv[grepl(time_0_contains, csv$`Biological sample category`)]$Sequence
    csv = csv[Sequence %in% seq_0]
  }
  annotation_table <- unique(csv[, .(Sequence, `Protein accession numbers`)])
  annotation_table <- na.omit(annotation_table)
  annotation_table$`Protein accession numbers` <- stringi::stri_split_fixed(annotation_table$`Protein accession numbers`, ",")
  annotation_table <- unique(annotation_table[, .(Proteins = unlist(`Protein accession numbers`)), by = Sequence])
  proteins <- unique(annotation_table$Proteins)
  require(Rcpp)
  sourceCpp("~MetaProfiler/fasta_reader.cpp")
  fasta <- unique(data.table::as.data.table(read_fasta(db, proteins)))
  fasta$sequence <- gsub("\n", "", fasta$sequence)
  annotation_table$`Protein Sequence` <- fasta[annotation_table$Proteins,sequence, on="accession"]
  annotation_table$`Protein Descriptions` <- fasta[annotation_table$Proteins, description, on="accession"]
  annotation_table <- na.omit(annotation_table, "Sequence")
  pro_seq = annotation_table$`Protein Sequence`
  pep_seq = annotation_table$Sequence
  pos = stringi::stri_locate_all_fixed(
    pro_seq, pep_seq
  )
  chars = strsplit(annotation_table$`Protein Sequence`, "")
  add = mapply(function(p, s) {
    s = c("[", s, "]")
    start = p[,1]
    end = p[,2]
    b = start
    a = end + 2
    aa_before = paste0(s[b], collapse = " ")
    aa_after = paste0(s[a], collapse = " ")
    start = paste0(start - 1, collapse = " ")
    end = paste0(end - 1, collapse = " ")
    list(aa_before = aa_before, aa_after = aa_after, start = start, end = end)
  }, pos, chars, SIMPLIFY = F)
  add = data.table::rbindlist(add)
  test = rowSums(add != "NA") != 0
  annotation_table[, colnames(add)] = add
  annotation_table = annotation_table[test]
  annotation_table = unique(annotation_table)
  col.names <- colnames(annotation_table)
  razor <- annotation_table[, c(lapply(.SD[, col.names[2:4], with = F], list), lapply(.SD[, col.names[5:8], with = F], paste0, collapse = " ")), by = Sequence]
  razor <- data.table::setnames(.accession_razor(razor, "Proteins", "Sequence", progress = T)[,1:2,with=F], c("Sequence", "Proteins"))
  annotation_table = annotation_table[razor, , on = c("Sequence", "Proteins")]
  csv[,col.names] <- annotation_table[csv$Sequence,,on="Sequence"]
  csv$`Observed m/z` <- as.numeric(gsub(",", "", csv$`Observed m/z`))
  csv[is.na(csv)] <- 0
  
  for (e in unique(csv$`Biological sample category`))
  {
    x = csv[`Biological sample category` == e]
    pro2id <- data.table::data.table(protein = unlist(x$Proteins), sequence = unlist(x$`Protein Sequence`),
                                     description = unlist(x$`Protein Descriptions`))
    pro2id <- unique(pro2id)
    pro2id$id <- 1:nrow(pro2id)
    map2pep <- rep(1:nrow(x), lengths(x$Proteins))
    ids <- pro2id[unlist(x$Proteins), id, on = "protein"]
    ids <- split(ids, map2pep)
    ids <- sapply(ids, function(x) paste0("PH_", x, collapse = " "))
    x$Protein_ids <- ids
    store(
      mz = x$`Observed m/z`,
      RT = as.numeric(stringi::stri_extract_first_regex(x$`Spectrum name`, "(?<=Elution: )[\\d\\.]+")) * 60,
      sequence = x$`Peptide sequence`,
      charge = x$`Spectrum charge`,
      aa_before = x$aa_before,
      aa_after = x$aa_after,
      start = x$start,
      end = x$end,
      document_id = e,
      Experiment = e,
      ids = ids,
      scan = x$`Spectrum name`,
      peptide_score = as.numeric(gsub("%", "", x$`Peptide identification probability`)),
      peptide_score_type = "Peptide identification probability",
      protein_score = as.numeric(gsub("%", "", x$`Protein identification probability`)),
      protein_score_type = "Protein identification probability",
      peptide_higher_score_better = T,
      protein_higher_score_better = T,
      pro2id = pro2id,
      output = output,
      db = db, 
      variable_modifications = "Oxidation (M)"
    )
  }
}

peptideProphet2idxml <- function(dir, db, output, pep_time_0_only = T)
{
  csv = data.table::fread(file.path(dir, "psm.tsv"))
  csv = csv[!grepl("DECOY_", Protein)]
  csv$`Modified Peptide` = gsub("\\[147\\]", "(Oxidation)", csv$`Modified Peptide`)
  csv$`Modified Peptide` = gsub("\\[43\\]", ".(Acetyl)", csv$`Modified Peptide`)
  csv$`Modified Peptide` = gsub("C", "C(Carbamidomethyl)", csv$`Modified Peptide`)
  csv$`Modified Peptide`[csv$`Modified Peptide` == ""] = gsub("C", "C(Carbamidomethyl)",
                                                              csv$Peptide)[csv$`Modified Peptide` == ""]
  p = csv$`Mapped Proteins` != ""
  csv$Protein[p] = paste0(csv$Protein[p], ", ", csv$`Mapped Proteins`[p])
  annotation_table <- unique(csv[, .(Sequence = Peptide, Protein)])
  annotation_table <- na.omit(annotation_table)
  annotation_table$Protein <- stringi::stri_split_fixed(annotation_table$Protein, ", ")
  annotation_table <- unique(annotation_table[, .(Protein = unique(unlist(Protein))), by = Sequence])
  
  protein <- unique(annotation_table$Protein)
  fasta <- unique(data.table::as.data.table(read_fasta(db, protein)))
  fasta$sequence <- gsub("\n", "", fasta$sequence)
  annotation_table$`Protein Sequence` <- fasta[annotation_table$Protein,sequence, on="accession"]
  annotation_table$`Protein Descriptions` <- fasta[annotation_table$Protein,description, on="accession"]
  annotation_table <- na.omit(annotation_table, "Sequence")
  pro_seq = annotation_table$`Protein Sequence`
  pep_seq = annotation_table$Sequence
  pos = stringi::stri_locate_all_fixed(
    pro_seq, pep_seq
  )
  chars = strsplit(annotation_table$`Protein Sequence`, "")
  add = mapply(function(p, s) {
    s = c("[", s, "]")
    start = p[,1]
    end = p[,2]
    b = start
    a = end + 2
    aa_before = paste0(s[b], collapse = " ")
    aa_after = paste0(s[a], collapse = " ")
    start = paste0(start - 1, collapse = " ")
    end = paste0(end - 1, collapse = " ")
    list(aa_before = aa_before, aa_after = aa_after, start = start, end = end)
  }, pos, chars, SIMPLIFY = F)
  add = data.table::rbindlist(add)
  test = rowSums(add != "NA") != 0
  annotation_table[, colnames(add)] = add
  annotation_table = annotation_table[test]
  annotation_table = unique(annotation_table)
  col.names <- colnames(annotation_table)
  annotation_table <- annotation_table[, c(lapply(.SD[, col.names[2:4], with = F], list), lapply(.SD[, col.names[5:8], with = F], paste0, collapse = " ")), by = Sequence]
  csv[,col.names] <- annotation_table[csv$Peptide,,on="Sequence"]
  f <- function(x)
  {
    if(is.character(x))
    {
      if(all(!grepl(",", x)))
      {
        return(x)
      }
      return(stringi::stri_split_regex(x, ",\\s*"))
    }
    return(x)
  }
  csv <- csv[, lapply(.SD, f)]
  csv$Experiment = gsub("\\.\\d+\\.\\d+\\.\\d+", "", csv$Spectrum)
  for (e in unique(csv$Experiment))
  {
    x = csv[Experiment == e]
    pro2id <- data.table::data.table(protein = unlist(x$Protein), sequence = unlist(x$`Protein Sequence`),
                                     description = unlist(x$`Protein Descriptions`))
    pro2id <- unique(pro2id)
    pro2id$id <- 1:nrow(pro2id)
    map2pep <- rep(1:nrow(x), lengths(x$Protein))
    ids <- pro2id[unlist(x$Protein), id, on = "protein"]
    ids <- split(ids, map2pep)
    ids <- sapply(ids, function(x) paste0("PH_", x, collapse = " "))
    x$Protein_ids <- ids
    x$Retention
    x$`Calculated M/Z`
    x$`Observed M/Z`
    store(
      mz = x$`Observed M/Z`,
      RT = x$Retention,
      sequence = x$`Modified Peptide`,
      charge = x$Charge,
      aa_before = x$aa_before,
      aa_after = x$aa_after,
      start = x$start,
      end = x$end,
      ids = ids,
      scan = x$Spectrum,
      peptide_score = x$`PeptideProphet Probability`,
      peptide_score_type = "PeptideProphet Probability",
      protein_score = 0,
      protein_score_type = "ProteinProphet Probability",
      peptide_higher_score_better = T,
      protein_higher_score_better = T,
      pro2id,
      output = output,
      db, 
      e, 
      e, 
      fixed_modifications = NULL, 
      variable_modifications = "Oxidation (M)"
    )
  }
}
