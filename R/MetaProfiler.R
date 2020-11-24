setClass(
  "MetaProfiler",
  slots = c(
    design = "data.table",
    data = "data.table",
    timepoints = "numeric",
    time_zero = "numeric",
    time_unit = "character",
    isotope = "character",
    incorporation_name = "character",
    incorporation_full_name = "character",
    incorporation_columns = "character",
    intensity_name = "character",
    intensity_full_name = "character",
    intensity_columns = "character",
    labeling_ratio_name = "character",
    labeling_ratio_full_name = "character",
    labeling_ratio_columns = "character",
    score_name = "character",
    score_full_name = "character",
    score_columns = "character",
    higher_score_better = "logical",
    as_percentage = "logical",
    peptide_column_no_PTMs = "character",
    peptide_column_PTMs = "character",
    accession_column = "character",
    pep2pro_peptide_column = "character",
    pep2pro_accession_column = "character",
    pep2func_peptide_column = "character",
    pep2func_function_columns = "character",
    pep2taxon_peptide_column = "character",
    rank_columns = "character",
    pro2func_accession_column = "character",
    pro2func_function_columns = "character",
    peptide_centric = "logical",
    annotate_with = "character",
    master_tbl = "data.table",
    annotation_table = "data.table",
    consts = "list",
    pdf = "list"
  ),
  prototype = list(
    design = data.table::data.table(),
    data = data.table::data.table(),
    timepoints=numeric(),
    time_zero = numeric(),
    time_unit = character(),
    isotope = character(),
    incorporation_name = character(),
    incorporation_columns = character(),
    intensity_name = character(),
    intensity_columns = character(),
    labeling_ratio_name = character(),
    labeling_ratio_columns = character(),
    score_name = character(),
    score_columns = character(),
    peptide_column_no_PTMs = character(),
    peptide_column_PTMs = character(),
    accession_column = character(),
    peptide_centric = logical(),
    annotate_with = character(),
    master_tbl = data.table::data.table(),
    annotation_table = data.table::data.table(),
    rank_columns = character(),
    pro2func_function_columns = character(),
    consts = list(),
    pdf = list()
  )
)

deconvolve = function(
  Object,
  ranges = rep(list(c(0,100)), 2),
  plot_condition = Object@design,
  nrow = NULL,
  ncol = NULL,
  filename = NULL,
  plot = TRUE,
  ...,
  distr = c(1,2,3)
) {
  
  plot_condition = unique(plot_condition[, sapply(plot_condition, function(col) {
    all(sapply(unique(col), function(val) {
      any(grepl(val, names(Object@pdf)))
    }))
  }), with = F])
  plot_condition = plot_condition[do.call(order, as.list(mget(names(plot_condition))))]
  if(!length(ncol) && length(nrow)) {
    ncol = ceiling(nrow(plot_condition)/nrow)
  }
  if(!length(nrow) && length(ncol)) {
    nrow = ceiling(nrow(plot_condition)/ncol)
  }
  if(!length(nrow) && !length(ncol)) {
    ncol = 3
    nrow = ceiling(nrow(plot_condition)/nrow)
  }
  obs = names(Object@pdf[[1]]$data)
  if(length(ranges) != length(obs)) {
    stop("Number of ranges given not equal to the number of variates in LFDR.")
  }
  
  if(!length(names(ranges))) {
    names(ranges) = obs
  } else if (!all(names(ranges) %in% obs)) {
    stop("Names provided for the ranges was not used as a variate for LFDR.")
  } else {
    ranges = ranges[obs]
  }
  
  den = function(n, i) {
    xs = lapply(ranges, function(r) {
      x = seq(r[1], r[2], length.out = n)
    })
    xs_ = lapply(xs, function(x) {
      x_ = x[-length(x)] + diff(x)/2
    })
    xs_c = mapply(function(x, d) {
      x_c = cut(d, x, include.lowest = T)
    }, xs, pdf[[i]]$data, SIMPLIFY = F)
    z1 = do.call(table, xs_c)
    denom = prod(sapply(lapply(xs, diff), '[[', 1))
    z1 = z1/(sum(z1)) * (1/denom)
    tbl = do.call(expand.grid, xs_)
    z2 = tryCatch(
      data.table::as.data.table(cbind(tbl, z = copula::dMvdc(as.matrix(tbl), pdf[[i]]$mvdc))),
      error = function(e) return(data.table::as.data.table(cbind(tbl, z = 0)))
    )
    # z2 = as.matrix(data.table::dcast(z2, x~y, value.var = "z")[,-1])
    z3 = data.table::as.data.table(cbind(tbl, z = ks::dkde(as.matrix(tbl), pdf[[i]]$kde0)))
    list(xs = xs, xs_ = xs_, xs_c = xs_c, z1 = z1, z2 = z2, z3 = z3)
  }
  
  to_opt = function(n, i) {
    penalty <- (round(n) - n)^2
    zs = den(n, i)
    # z3 = as.matrix(data.table::dcast(z3, x~y, value.var = "z")[,-1])
    sum(abs(as.data.frame(zs$z1)$Freq * 100 - (zs$z2$z  * (1-pdf[[i]]$prior_prob) + zs$z3$z * pdf[[i]]$prior_prob) * 100)^2)
  }
  pdf = Object@pdf[apply(plot_condition, 1, function(x) {
    y = paste0("\\b", names(plot_condition), " ", x, "\\b")
    which(rowSums(sapply(y, grepl, x = names(Object@pdf))) == length(x))
  })]
  
  grobs = NULL
  
  for(i in seq_len(length(pdf))) {
    from = 10
    to = 50
    n = from:to
    n = n[which.min(sapply(n, to_opt, i = i))]
    zs = den(n, i)
    top = NULL
    if(length(obs) == 2) {
      top = ggplotify::as.grob(function() {
        zs2 = den(512, i)
        true1 = zs2$z2[, .(z = sum(z)), by = c(obs[[1]])]
        true1$z = true1$z/sum(true1$z) * 1/diff(zs2$xs[[c(obs[[1]])]])
        false1 = zs2$z3[, .(z = sum(z)), by = c(obs[[1]])]
        false1$z = false1$z/sum(false1$z) * 1/diff(zs2$xs[[c(obs[[1]])]])
        
        true2 = zs2$z2[, .(z = sum(z)), by = c(obs[[2]])]
        true2$z = true2$z/sum(true2$z) * 1/diff(zs2$xs[[c(obs[[2]])]])
        false2 = zs2$z3[, .(z = sum(z)), by = c(obs[[2]])]
        false2$z = false2$z/sum(false2$z) * 1/diff(zs2$xs[[c(obs[[2]])]])
        panelfirst <- function(pmat) {
          rs1 = rescale(c(false1$z*(pdf[[i]]$prior_prob), true1$z*(1-pdf[[i]]$prior_prob),
                        false1$z*(pdf[[i]]$prior_prob) + true1$z*(1-pdf[[i]]$prior_prob)), 
                        to = range(zs$z1))
          rs1 = split(rs1, rep(c("false", "true", "joint"), each = length(true1$z)))
          XY <- plot3D::trans3D(y = rep(min(true2[[1]]) - max(true2[[1]])*0.5, 511), x = true1[[1]],
                        z = rs1$true, pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "#1E88E5", type = "l", lwd = 1, add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(y = rep(min(false2[[1]]) - max(false2[[1]])*0.5, 511), x = false1[[1]],
                                z = rs1$false, pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "#D81B60", type = "l", lwd = 1, add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(y = rep(min(false2[[1]]) - max(false2[[1]])*0.5, 511), x = false1[[1]],
                                z = rs1$joint,
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, add = TRUE, colkey = FALSE)
          
          
          
          rs2 = rescale(c(false2$z*(pdf[[i]]$prior_prob), true2$z*(1-pdf[[i]]$prior_prob),
                          false2$z*(pdf[[i]]$prior_prob) + true2$z*(1-pdf[[i]]$prior_prob)), 
                        to = range(zs$z1))
          rs2 = split(rs2, rep(c("false", "true", "joint"), each = length(true2$z)))
          XY <- plot3D::trans3D(y = true2[[1]], x = rep(min(true1[[1]]) - max(true2[[1]])*0.5, 511),
                        z = rs2$true, pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "#1E88E5", type = "l", lwd = 1, add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(y = false2[[1]], x = rep(min(false1[[1]]) - max(false2[[1]])*0.5, 511),
                                z = rs2$false, pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "#D81B60", type = "l", lwd = 1, add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(y = false2[[1]], x = rep(min(false1[[1]]) - max(false2[[1]])*0.5, 511),
                                z = rs2$joint, pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, add = TRUE, colkey = FALSE)
          
          
          XY <- plot3D::trans3D(y = c(min(false2[[1]]),min(false2[[1]])),
                                x = c(min(false2[[1]]) - max(false2[[1]])*0.5, min(false2[[1]])),
                                z = c(0,0),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3,
                            add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(y = c(max(false2[[1]]),max(false2[[1]])),
                                x = c(min(false2[[1]]) - max(false2[[1]])*0.5, min(false2[[1]])),
                                z = c(0,0),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3,
                            add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(x = c(max(false2[[1]]),max(false2[[1]])),
                                y = c(min(false2[[1]]) - max(false2[[1]])*0.5, min(false2[[1]])),
                                z = c(0,0),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3,
                            add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(x = c(min(false2[[1]]),min(false2[[1]])),
                                y = c(min(false2[[1]]) - max(false2[[1]])*0.5, min(false2[[1]])),
                                z = c(0,0),
                                pmat = pmat)
          
          
          XY <- plot3D::trans3D(y = c(min(false2[[1]]),min(false2[[1]])),
                                x = c(min(false2[[1]]) - max(false2[[1]])*0.5, min(false2[[1]])),
                                z = c(max(zs$z1),max(zs$z1)),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3,
                            add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(y = c(max(false2[[1]]),max(false2[[1]])),
                                x = c(min(false2[[1]]) - max(false2[[1]])*0.5, min(false2[[1]])),
                                z = c(max(zs$z1),max(zs$z1)),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3,
                            add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(x = c(max(false2[[1]]),max(false2[[1]])),
                                y = c(min(false2[[1]]) - max(false2[[1]])*0.5, min(false2[[1]])),
                                z = c(max(zs$z1),max(zs$z1)),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3,
                            add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(x = c(min(false2[[1]]),min(false2[[1]])),
                                y = c(min(false2[[1]]) - max(false2[[1]])*0.5, min(false2[[1]])),
                                z = c(max(zs$z1),max(zs$z1)),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3,
                            add = TRUE, colkey = FALSE)
          

          XY <- plot3D::trans3D(y = range(false2[[1]]), x = rep(min(false1[[1]]) - max(false2[[1]])*0.5, 2),
                                z = c(0,0),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3, 
                            add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(y = range(false2[[1]]), x = rep(min(false1[[1]]) - max(false2[[1]])*0.5, 2),
                                z = c(max(zs$z1),max(zs$z1)),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3, 
                            add = TRUE, colkey = FALSE)
          
          XY <- plot3D::trans3D(y = rep(min(false2[[1]]) - max(false2[[1]])*0.5, 2), x = range(false1[[1]]),
                                z = c(0,0),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3, 
                            add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(y = rep(min(false2[[1]]) - max(false2[[1]])*0.5, 2), x = range(false1[[1]]),
                                z = c(max(zs$z1),max(zs$z1)),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3, 
                            add = TRUE, colkey = FALSE)
          
          XY <- plot3D::trans3D(y = rep(max(false2[[1]]), 2),
                                x = rep(min(false1[[1]]) - max(false1[[1]])*0.5, 2),
                                z = c(0,max(zs$z1)),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3, 
                            add = TRUE, colkey = FALSE)
          
          XY <- plot3D::trans3D(y = rep(min(false2[[1]]), 2),
                                x = rep(min(false1[[1]]) - max(false1[[1]])*0.5, 2),
                                z = c(0,max(zs$z1)),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3, 
                            add = TRUE, colkey = FALSE)
          
          XY <- plot3D::trans3D(x = rep(min(false1[[1]]), 2),
                                y = rep(min(false2[[1]]) - max(false2[[1]])*0.5, 2),
                                z = c(0,max(zs$z1)),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3, 
                            add = TRUE, colkey = FALSE)
          
          XY <- plot3D::trans3D(x = rep(max(false1[[1]]), 2),
                                y = rep(min(false2[[1]]) - max(false2[[1]])*0.5, 2),
                                z = c(0,max(zs$z1)),
                                pmat = pmat)
          plot3D::scatter2D(XY$x, XY$y, col = "black", type = "l", lwd = 1, lty = 3, 
                            add = TRUE, colkey = FALSE)
          XY <- plot3D::trans3D(x = max(false1[[1]]),
                                y = min(false2[[1]]) - max(false2[[1]])*0.6,
                                z = diff(c(0,max(zs$z1)))/2,
                                pmat = pmat)
        }
        plot3D::hist3D(
          x = zs$xs_[[1]], y = zs$xs_[[2]], z = zs$z1, alpha = 0.5, prob = T, colkey = F, border = "black",
          phi = 20, theta = 135, col = matrix("white", ncol = ncol(zs$z1) - 1, nrow = nrow(zs$z1) - 1),
          colvar = NULL,  xlab = Object@incorporation_name, ylab = Object@labeling_ratio_name,
          zlab = "density", main = names(pdf)[i], panel.first = panelfirst, expand = 1, r = 10
        )
        z = zs$z3
        z$z = zs$z3$z*pdf[[i]]$prior_prob + zs$z2$z  * (1 - pdf[[i]]$prior_prob)
        z = as.matrix(data.table::dcast(z, get(obs[1])~get(obs[2]), value.var = "z")[,-1])
        plot3D::persp3D(x = zs$xs_[[1]], y = zs$xs_[[2]], z = z, curtain = T,  border = "purple", 
                        add = T, colkey = F, alpha = 0.3, main = NULL, pmat)
        
      })
    }
    grobs = c(grobs, list(top))
    writeLines(paste0(i,"/",length(pdf)))
  }
  gg = gridExtra::arrangeGrob(grobs = grobs, ncol = ncol, nrow = nrow)
  save_plot(filename = filename, plot = gg, ...)
}
get_cols = function(Object, var, x = Object@timepoints, time_unit = Object@time_unit, return_tbl = F) {
  if(!is.list(x))
    x = list(x)
  if(length(x) == 1)
    x = rep(x, length(var))
  col = NULL
  for(i in seq_along(var)) {
    nt = length(x[[i]])
    nv = length(var[[i]])
    col = c(col, paste0(rep(var[[i]], each = nt), " (", time_unit, " ", rep(x[[i]], nv), ")"))
  }
  if(return_tbl)
    return(Object@master_tbl[, ..col])
  col
}


model = function(Object = NULL, x, var, take_mean = FALSE, subset = NULL, consts = NULL) {
  if(length(Object)) {
    consts = Object@consts[[var]]
  } else if (!length(consts)) {
    stop("Object is missing and no consts table is provided")
  }
  if(!length(consts)) {
    stop(paste0("Curve fitting was not performed on variable `", var, "`."))
  }
  eq <- unique(consts$equation)
  if(eq == 4) {
    consts <- consts
    value <- lapply(consts$k, predict, x)
    value <- as.data.table(do.call(rbind, lapply(consts$k, predict, x))) 
    value$Peptides <- rep(consts[, Peptides], each = length(x))
    value <- as.matrix(dcast(value, Peptides ~ z, value.var = "fit"), rownames = TRUE)
  } else {
    if(!length(subset)) {
      subset <- 1:nrow(consts)
    }
    if(length(x) > 1 && !is.matrix(x))
    {
      x <- matrix(x, nrow = nrow(consts[subset,]), ncol = length(x), byrow = TRUE)
    }
    if (eq == 1) {
      k = consts[[1]][subset]
      conc = consts[[2]][subset]
      value = (1 - exp(-k * x)) * conc
      is_percentage = .check_if_percantage_or_fraction(value)
      if(length(Object) && !Object@as_percentage && is_percentage) {
        value = value / 100
      }
      if(length(Object) && Object@as_percentage && !is_percentage) {
        value = value * 100
      }
      value
    }
    else if(eq == 2) {
      kbi = consts[[1]][subset]
      k0a = consts[[2]][subset]
      conc = consts[[3]][subset]
      v = kbi
      yv = kbi/(k0a - kbi)
      ykbi = k0a/(kbi - k0a)
      value = (1 + yv * exp(-v * x) + ykbi * exp(-kbi * x)) * conc
      is_percentage = .check_if_percantage_or_fraction(value)
      if(length(Object) && !Object@as_percentage && is_percentage) {
        value = value / 100
      }
      if(length(Object) && Object@as_percentage && !is_percentage) {
        value = value * 100
      }
      value
    } else if (eq == 3) {
      k0a = consts[[1]][subset]
      kbi = consts[[2]][subset]
      kst = consts[[3]][subset]
      kbt = consts[[4]][subset]
      conc = consts[[5]][subset]
      u = ((kst + k0a + kbt) - sqrt((kst + k0a + kbt)^2 - (4 * k0a*kbt))) / 2
      v = ((kst + k0a + kbt) + sqrt((kst + k0a + kbt)^2 - (4 * k0a*kbt))) / 2
      yu = k0a * kbi*(u - kbt) / ((u - v)*(u - kbi)*u)
      yv = k0a * kbi*(v - kbt) / ((v - u)*(v - kbi)*v)
      ykbi = k0a * (kbi - kbt) / ((u - kbi)*(v - kbi))
      value = (1 + yu * exp(-u * x) + yv * exp(-v * x) + ykbi * exp(-kbi * x))*conc
      is_percentage = .check_if_percantage_or_fraction(value)
      if(length(Object) && !Object@as_percentage && is_percentage) {
        value = value / 100
      }
      if(length(Object) && Object@as_percentage && !is_percentage) {
        value = value * 100
      }
      value
    }
  }
  if(take_mean) {
    value = colMeans(as.matrix(value))
  }
  return(value)
}

#' @title Create MetaProfiler class object from protein-SIP results.
#' @description Converts the results obtained from protein-SIP experiment into a MetaProfiler class. If using result 
#' files, then the extension must be either be a csv or tsv. 
#' @param design A \code{data.frame} or \code{data.table} containing the experimental design. Each row must 
#' correspond to the file paths listed in \code{data} or in the file column of table \code{design}. Additionally, 
#' a timepoint column with the same name as \code{time_unit} must be present.
#' @param data Either a \code{data.frame}/\code{data.table} or a character vector of the list of files containing 
#' the result from the protein-SIP experiment. Must be tab or comma delimited. Can be set to \code{NULL} 
#' if table \code{design} contains a column with the filepaths of the result files.
#' @param time_unit The unit for the timepoints.
#' @param time_zero Numeric value denoting the timepoint when the diet was switched. Defaults to zero.
#' @param isotope A character value specifying the stable isotope. Should correspond to one of the elements in the 
#' periodic table.  
#' @param peptide_centric Logical value specifying if analysis is done at the peptide or protein level.
#' @param incorporation_name A character value denoting the name of the incorporation value. 
#' If \code{incorporation_columns} is set to \code{"auto"}, then the function will look for columns 
#' that start with the character value followed by a unique identifier (i.e. "RIA 1", "RIA 2", "RIA 3", ...). 
#' See details for the difference between incorporation and labeling ratio.
#' @param incorporation_columns A character vector detailing the names of the columns containing the incorporation 
#' values. Can be set to \code{"auto"} if \code{incorporation_name} is specified. Can also be set to \code{NULL} 
#' if no incorporation value was measured.
#' @param intensity_name Silmilar to \code{incorporation_name}, but with the intensity values instead.
#' @param intensity_columns Silmilar to \code{incorporation_columns}, but with the intensity columns instead. 
#' Can be set to \code{NULL} if labeling ratio was measured instead.
#' @param labeling_ratio_name Silmilar to \code{incorporation_name}, but with the labeling ratio values instead. 
#' See details for the difference between incorporation and labeling ratio.
#' @param labeling_ratio_columns Silmilar to \code{incorporation_columns}, but with the labeling ratio columns 
#' instead. If set to \code{NULL}, then labeling ratio is calculated from intensity.
#' @param labeling_ratio_threshold A numeric value that specifies the minimum labeling ratio needed for the 
#' heavy peptide or protein to be kept for downstream analysis.
#' @param score_name Silmilar to \code{incorporation_name}, but with the heavy peptide identification score values 
#' instead.
#' @param score_columns Silmilar to \code{incorporation_columns}, but with the heavy peptide identification score 
#' columns instead. Can be set to \code{NULL} if score was not measured.
#' @param score_threshold A numeric value that specifies the minimum or maximum score needed for the heavy peptide or 
#' protein to be kept for downstream analysis.
#' @param peptide_column_PTMs A character value. Specifies the name of the column containing the peptide sequence 
#' with post translational modifications (PTMs). If set to \code{"guess"}, then the function will guess the column 
#' based off the headers and the characters contained in the column.
#' @param peptide_column_no_PTMs Similar to \code{peptide_column_PTMs}, but instead with peptide sequences without 
#' post translational modifications (PTMs). If set to \code{"guess"}, but the function only detects the column 
#' containing PTMs, then the function will add a new column containing the sequences from \code{peptide_column_PTMs}, 
#' but with modications removed. This also aplies to when the value is set to \code{NULL}. When removing PTMs, 
#' the function assumes that the name of the modications in the sequences follow UniProt convention.
#' @param accession_column A character value which specifies the name of the column containing the protein accession 
#' IDs. Can be set to \code{NULL} if \code{pep2pro} is provided.
#' @param accession_pattern A string regex. Only used when \code{accession_column} is set to \code{"guess"}. In this 
#' case, the function will find columns containing regex strings from \code{accession_pattern}. 
#' Defaults to standard naming convention patterns for UniProt IDs.
#' @param compute_razor_protein If set to \code{TRUE}, then the razor protein for each peptide will be computed.
#' See details for information about razor proteins.
#' @param pep2pro A character value for the filename or a \code{data.frame}/\code{data.table} with the peptide to protein table.
#' Optional if \code{accession_column} is provided.
#' @param pep2pro_peptide_column A character value for the peptide column in \code{pep2pro}. 
#' By default, the function will guess the columns similarly to \code{peptide_column_PTMs} or 
#' \code{peptide_column_no_PTMs}, depending on \code{annotate_by_peptide_with_PTMs}.
#' @param pep2pro_accession_column A character value for the protein accessions IDs column in \code{pep2pro}. 
#' By default, the function will guess the column similarly to how \code{accession_column} is guessed.
#' @param pep2pro_accession_pattern See \code{accession_pattern}.
#' @param pep2taxon A character value for the filename or a \code{data.frame}/\code{data.table} with the peptide to taxon table.
#' @param pep2taxon_peptide_column A character value for the peptide column in \code{pep2taxon}. 
#' By default, the function will guess the columns similarly to \code{peptide_column_PTMs} or 
#' \code{peptide_column_no_PTMs}, depending on \code{annotate_by_peptide_with_PTMs}.
#' @param rank_columns A character vector for the phylogenetic rank columns in \code{pep2taxon}. By default,
#' the function will look for columns with the names `superkingdom`, `kingdom`, `subkingdom`, `superphylum`, 
#' `phylum`, `subphylum`, `superclass`, `class`, `subclass`, `infraclass`, `superorder`, `order`, `suborder`, 
#' `infraorder`, `parvorder`, `superfamily`, `family`, `subfamily`, `tribe`, `subtribe`, `genus`, `subgenus`, 
#' `species group`, `species subgroup`, `species`, `subspecies`, `varietas`, and `forma`.
#' @param pro2func A character value for the filename or a \code{data.frame}/\code{data.table} with the protein to function table.
#' @param pro2func_accession_column A character value for the protein accessions IDs column in \code{pro2func}. 
#' By default, the function will guess the column similarly to how \code{accession_column} is guessed.
#' @param pro2func_accession_pattern See \code{accession_pattern}.
#' @param pro2func_function_columns A character vector for the protein function columns in \code{pro2func}. 
#' By default, the function will guess the columns by looking for columns starting with common functional annotation databases such as 
#' KEGG, BRITE, GO, COG, and NOG.
#' @param annotate_by_peptide_with_PTMs A logical value that specifies whether functional and taxonomic annotation 
#' is done using peptide sequences with PTMs, \code{TRUE}, or without, \code{FALSE}. Only used when peptide_centric 
#' is set to \code{TRUE}.
#' @param feature_type_column A character value specifying the name of the column containing the type of feature 
#' the heavy peptide was identified/quantified from. If using results from MetaProSIP, the reference used for heavy 
#' peptide identification can be either from a feature (i.e. the group of peaks in the retention time and mass to 
#' charge ratio dimension belonging to a single peptide entity) or from a pseudo-feature (i.e. the theoretical position 
#' of the unlabeled feature using sequence information only).
#' @param feature_type A character vector for the types of feature contained in \code{feature_type_column}. 
#' Defaults to \code{"feature"} and \code{"id"} (i.e. the pseudo-feature).
#' @param as_percentage Should incorporation and labeling ratio values be presented as percentages?
#' @param progress If \code{TRUE}, then progress is printed.
#' @param trace If \code{TRUE}, then the log is printed.
#' @details 
#' # Labeling Ratio
#' \code{labeling_ratio_name} denotes the relative ratio between the light peptide and the estimated 
#' intensity of the heavy peptide. When using an unlabeled protein-spike in, LR measures the proportion of proteins 
#' that are produced using the heavy stable isotope relative to the protein at \code{time_zero}. By taking this 
#' measure over time, the rate of newly synthesized proteins that incorporate the stable isotope can be estimated. 
#' Labeling ratio is calculated using equation:
#' \deqn{(Ih)/(Il+Ih),}
#' where \eqn{Il} is the sum of the intensities of the unlabeled peptides or proteins and \eqn{Ih} is the sum of the 
#' intensities of the heavy peptides or proteins.
#' # Elemental Flux
#' \code{incorporation_name} describes the elemental flux of the isotope, which is measured using the average proportion of the stable isotope 
#' incorporated in a peptide of interest. By characterizing the functional and taxonomic origin of the peptide, 
#' it gives insight on how and where the stable isotopic substrate is being converted into biomass. Thus, measuring 
#' the elemental flux can predict where this substrate is limited. Incorporation is calculated using equation:
#' \deqn{(H-M)/(FALSE-M),}
#' where \eqn{H} is the m/z position at the center of the predicted isotopic distribution of the heavy peptide, 
#' \eqn{M} is the monoisotopic peak of the light peptide, and \eqn{FALSE} is the m/z position of the fully labeled peptide.
#' @return Returns an object of class MetaProfiler.
#' @examples
#' 
MetaProfiler <- function(design,
                         data,
                         time_unit,
                         time_zero = 0,
                         isotope = "N",
                         peptide_centric = TRUE,
                         light_peptide = TRUE,
                         as_percentage = TRUE,
                         
                         incorporation_name = "RIA",
                         incorporation_full_name = "Relative Isotopic Abundance",
                         incorporation_columns = "auto",
                         
                         intensity_name = "INT",
                         intensity_full_name = "Intensity",
                         intensity_columns = "auto",
                         
                         labeling_ratio_name = "LR",
                         labeling_ratio_full_name = "Labeling Ratio",
                         labeling_ratio_columns = NULL,
                         labeling_ratio_threshold = NA,
                         
                         score_name = "Cor",
                         score_full_name = "Correlation",
                         higher_score_better = TRUE,
                         score_columns = "auto", 
                         score_threshold = NA,
                         
                         peptide_column_PTMs = "guess", 
                         peptide_column_no_PTMs = "guess",
                         accession_column = "guess",
                         accession_pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",
                         compute_razor_protein = FALSE,
                         
                         pep2pro = NULL,
                         pep2pro_peptide_column = "guess",
                         pep2pro_accession_column = "guess",
                         pep2pro_accession_pattern = accession_pattern,
                         
                         pep2func = NULL,
                         pep2func_peptide_column = "guess",
                         pep2func_function_columns = "guess",
                         
                         pep2taxon = NULL,
                         pep2taxon_peptide_column = "guess",
                         rank_columns = "guess",
                         main_ranks_only = T,
                         
                         pro2func = NULL,
                         pro2func_accession_column = "guess",
                         pro2func_accession_pattern = accession_pattern,
                         pro2func_function_columns = "guess",
                         annotate_by_peptide_with_PTMs = FALSE,
                         
                         feature_type_column = NULL,
                         feature_type = c("feature", "id"),
                         
                         progress = TRUE, 
                         trace = TRUE,
                         raw = F,
                         ...) {
  Object <- new("MetaProfiler")
  Object@design <- data.table::as.data.table(design)
  Object@time_unit <- time_unit
  Object@time_zero <- time_zero
  Object@higher_score_better = higher_score_better
  Object@peptide_centric = peptide_centric
  Object@isotope = isotope
  Object@as_percentage = as_percentage
  if(length(incorporation_name))
    Object@incorporation_name = incorporation_name
  if(length(incorporation_full_name))
    Object@incorporation_full_name = incorporation_full_name
  if(length(intensity_name))
    Object@intensity_name = intensity_name
  if(length(intensity_full_name))
    Object@intensity_full_name = intensity_full_name
  if(length(labeling_ratio_name))
    Object@labeling_ratio_name = labeling_ratio_name
  if(length(labeling_ratio_full_name))
    Object@labeling_ratio_full_name = labeling_ratio_full_name
  if(length(score_name))
    Object@score_name = score_name
  if(length(score_full_name))
    Object@score_full_name = score_full_name
  if(length(incorporation_columns) && (incorporation_columns != "auto"))
    Object@incorporation_columns = incorporation_columns
  if(length(intensity_columns) && (intensity_columns != "auto"))
    Object@intensity_columns = intensity_columns
  if(length(labeling_ratio_columns) && (labeling_ratio_columns != "auto"))
    Object@labeling_ratio_columns = labeling_ratio_columns
  if(length(score_columns) && (score_columns != "auto"))
    Object@score_columns = score_columns
  
  if(!any(colnames(Object@design) == Object@time_unit)) {
    stop(paste0("No column named `", Object@time_unit, "` in the table `design`. Please add/rename the column for the timepoints or change variable `time_unit` so that it matches the timepoint column in table `design`."))
  }
  if(!is.numeric(Object@design[[Object@time_unit]])) {
    Object@design[[Object@time_unit]] <- as.numeric(Object@design[[Object@time_unit]])
    if(all(is.na(Object@design[[Object@time_unit]]))) {
      stop(paste0("Could not convert column `", Object@time_unit, "` to numeric. Are there any non-numeric characters in the column?"))
    }
  }
  check = unlist(Object@design[, lapply(.SD, function(x) {
    test <- try(all(file.exists(x)), silent = TRUE)
    ifelse(inherits(test, "try-error"), FALSE, test)
  })])
  if(missing(data)) {
    data <- NULL
  }
  if(!inherits(data, "data.frame")) {
    if(length(data) && is.character(data)) {
      if (nrow(Object@design) != length(data)) {
        stop("Make sure to provide the experimental conditions (e.g. sample names, timepoints) for all the files listed in `data`")
      } else {
        if(length(Object@design) > 1) {
          nm <- names(check)[!check]
          Object@design <- Object@design[, ..nm]
        }
        Object@design <- cbind(filenames = data, Object@design)
      }
    } else if (sum(check) == 1) {
      idx <- order(check, decreasing = TRUE)
      Object@design <- Object@design[,idx,with=FALSE]
    } else if (any(check)) {
      stop("More than one column in table `design` contains filepaths. Initialize `data` with the appropriate filepaths instead.")
    }
    readr_params <- list(...)
    options("readr.show_progress" = FALSE)
    f <- function(...)
    {
      args <- list(...)
      if(!file.exists(args[[1]])) {
        stop("Either the first column in `design` does not contain the filenames or that file \"", x, "\" does not exist.")
      }
      if(!any(names(readr_params) == "delim")) {
        readr_params$delim = "\t"
      }
      if(progress) {
        cat(paste0("reading ", gsub("[\\S\\s]+\\/", "", args[[1]], perl = TRUE), " (", paste0(colnames(Object@design)[-1], " = ", unlist(args)[-1], collapse = ", "),")...\n"))
      }
      shush(data <- do.call(readr::read_delim, c(list(file = args[[1]]), readr_params)))
      if(progress) {
        cat("done.\n")
      }
      data <- data.table::as.data.table(setNames(cbind(args[-1], data), c(colnames(Object@design)[-1], colnames(data))))
      data
    }
    args <- c(list(FUN = f), as.list(Object@design), list(SIMPLIFY = FALSE))
    data <- data.table::rbindlist(do.call(mapply, args))
  }
  headers = colnames(data)
  if(length(incorporation_columns) == 1 && incorporation_columns == "auto") {
    incorporation_columns = .check_columns(incorporation_name, headers)
  }
  if(length(intensity_columns) == 1 && intensity_columns == "auto") {
    intensity_columns = .check_columns(intensity_name, headers)
  }
  if(length(score_columns) == 1 && score_columns == "auto") {
    score_columns = .check_columns(score_name, headers)
  }
  if(length(labeling_ratio_columns) == 1 && labeling_ratio_columns == "auto") {
    labeling_ratio_columns = .check_columns(labeling_ratio_name, headers)
  }
  
  if(length(incorporation_columns))
    Object@incorporation_columns = incorporation_columns
  if(length(intensity_columns))
    Object@intensity_columns = intensity_columns
  if(length(score_columns))
    Object@score_columns = score_columns
  if(length(labeling_ratio_columns))
    Object@labeling_ratio_columns = labeling_ratio_columns
  group = c(incorporation_columns, intensity_columns, score_columns, labeling_ratio_columns)
  if (length(group)) {
    data[,group] <- data[, lapply(.SD, as.numeric), .SD = group]
    if (length(incorporation_columns) || length(labeling_ratio_columns)) {
      data[, c(incorporation_columns, labeling_ratio_columns)] = data[, lapply(.SD, function(x) {
        is_percentage <- .check_if_percantage_or_fraction(x)
        if(!length(is_percentage)) {
          return(as.numeric(x))
        } else if (is_percentage) {
          if(as_percentage) {
            return(x)
          } else {
            return(x * 100)
          }
        } else {
          if(as_percentage) {
            return(x * 100)
          } else {
            return(x)
          }
        }
      }), .SDcols = c(incorporation_columns, labeling_ratio_columns)]
    }
  }
  if(length(Object@design) > 1) {
    Object@design <- Object@design[, !check, with = FALSE]
  }
  
  
  res <- .guess_columns(
    data,
    peptide_column_PTMs,
    peptide_column_no_PTMs,
    accession_column,
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
  )
  
  Object@annotate_with = res$annotate_with
  Object@data = res$data
  if(length(res$data_peptide_column_PTMs))
    Object@peptide_column_PTMs = res$data_peptide_column_PTMs
  if(length(res$data_peptide_column_no_PTMs))
    Object@peptide_column_no_PTMs = res$data_peptide_column_no_PTMs
  if(length(res$data_accession_column))
    Object@accession_column = res$data_accession_column
  
  if(length(res$pep2pro_peptide_column))
    Object@pep2pro_peptide_column = res$pep2pro_peptide_column
  if(length(res$pep2pro_accession_column))
    Object@pep2pro_accession_column = res$pep2pro_accession_column
  
  if(length(res$pep2func_peptide_column))
    Object@pep2func_peptide_column = res$pep2func_peptide_column
  if(length(res$pep2func_function_columns))
    Object@pep2func_function_columns = res$pep2func_function_columns
  
  if(length(res$pep2taxon_peptide_column))
    Object@pep2taxon_peptide_column = res$pep2taxon_peptide_column
  if(length(res$rank_columns))
    Object@rank_columns = res$rank_columns
  
  if(length(res$pro2func_accession_column))
    Object@pro2func_accession_column = res$pro2func_accession_column
  if(length(res$pro2func_function_columns))
    Object@pro2func_function_columns = res$pro2func_function_columns
  
  if(length(feature_type_column)) {
    Object@data <- Object@data[feature_type, , on = feature_type_column]
  }
  Object@timepoints = sort(unique(Object@data[,get(Object@time_unit)]))
  Object@timepoints = Object@timepoints[Object@timepoints != Object@time_zero]
  Object <- .make_annotation_table(
    Object = Object,
    compute_razor_protein = compute_razor_protein,
    pep2pro = res$pep2pro,
    pep2func = res$pep2func,
    pep2taxon = res$pep2taxon,
    pro2func = res$pro2func,
    trace = trace
  )
  if(!raw) {
    Object <- .get_data(Object = Object, light_peptide = light_peptide, trace = trace, score_threshold = score_threshold, labeling_ratio_threshold = labeling_ratio_threshold)
  }
  options("readr.show_progress" = TRUE)
  Object
}

.get_data <- function(Object, light_peptide, trace, score_threshold, labeling_ratio_threshold) {
  Object@data = data.table::copy(Object@data)
  if (length(Object@incorporation_columns) && (length(Object@intensity_columns) || length(Object@labeling_ratio_columns))) {
    light_cols <- na.omit(c(Object@incorporation_columns[1], Object@intensity_columns[1], Object@score_columns[1], Object@labeling_ratio_columns[1]))
    heavy_cols <- list(Object@incorporation_columns[-1], Object@intensity_columns[-1], Object@score_columns[-1], Object@labeling_ratio_columns[-1])
    keep = !sapply(heavy_cols, function(x) {
      length(x) == 0
    })
    heavy_cols = heavy_cols[keep]
    value.name = unlist(list(Object@incorporation_name, Object@intensity_name, Object@score_name, Object@labeling_ratio_name)[keep])
    seq <- unique(Object@data[[Object@annotate_with]])
    seq <- setNames(strsplit(seq, "", fixed = TRUE), seq)
    seq <- lapply(seq, function(x) setNames(as.list(x), as.character(seq_along(x))))
    seq <- data.table::rbindlist(seq, fill = TRUE, idcol = "pep")
    seq[,colnames(seq)[-1]] <- seq[,lapply(.SD[,-1], function(x) AA_tbl[x, , on = "AA"][[Object@isotope]])]
    seq <- seq[, sum(unlist(.SD), na.rm = TRUE), by = pep]
    Object@data$isotopes <- seq[Object@data[[Object@annotate_with]], V1, on = "pep"]
    incorporation_value <- Object@data[,Object@incorporation_columns,with=FALSE]
    is_percentage = .check_if_percantage_or_fraction(incorporation_value)
    if(is_percentage) {
      incorporation_value <- incorporation_value/100
    }
    test <- rowSums(incorporation_value > (1/Object@data$isotopes), na.rm = TRUE) > 0
    if(!light_peptide) {
      swap <- is.na(incorporation_value[[2]]) & test
      incorporation_value[[2]][swap] = incorporation_value[[1]][swap]
      incorporation_value[[1]][swap] = 0
      if(is_percentage)
        Object@data[,Object@incorporation_columns] <- incorporation_value * 100
      else
        Object@data[,Object@incorporation_columns] <- incorporation_value
      if(length(Object@intensity_columns)){
        Object@data[[Object@intensity_columns[2]]][swap] = Object@data[[Object@intensity_columns[1]]][swap]
        Object@data[[Object@intensity_columns[1]]][swap] = 0
      }
      if(length(Object@score_columns)){
        Object@data[[Object@score_columns[2]]][swap] = Object@data[[Object@score_columns[1]]][swap]
        Object@data[[Object@score_columns[1]]][swap] = 0
      }
      if(length(Object@labeling_ratio_columns)){
        Object@data[[Object@labeling_ratio_columns[2]]][swap] = Object@data[[Object@labeling_ratio_columns[1]]][swap]
        Object@data[[Object@labeling_ratio_columns[1]]][swap] = 0
      }
      if(length(Object@score_columns)){
        Object@data[[Object@score_columns[2]]][swap] = Object@data[[Object@score_columns[1]]][swap]
        if(Object@higher_score_better) {
          Object@data[[Object@score_columns[1]]][swap] = max(Object@data[,Object@score_columns,with=FALSE],na.rm = TRUE)
        } else {
          Object@data[[Object@score_columns[1]]][swap] = min(Object@data[,Object@score_columns,with=FALSE],na.rm = TRUE)
        }
      }
    }
    test <- test & (rowSums(incorporation_value[,1] <= (1/Object@data$isotopes), na.rm = TRUE) > 0)
    Object@data <- Object@data[test]
    measurements <- c(Object@incorporation_columns, Object@intensity_columns, Object@score_columns, Object@labeling_ratio_columns)
    group <- colnames(Object@data)[!(colnames(Object@data) %in% c(Object@incorporation_columns, Object@intensity_columns, Object@score_columns, Object@labeling_ratio_columns))]
    col <-  colnames(Object@data)[colnames(Object@data) %in% c(group, light_cols)]
    Object@data[, measurements] <- Object@data[, lapply(.SD, as.numeric), .SDcols = measurements]
    Object@data <- data.table::melt(Object@data, id.vars = col, measure = heavy_cols, value.name = value.name, na.rm = TRUE)[,-"variable"]
    if(length(Object@labeling_ratio_columns) == 1) {
      if(length(Object@labeling_ratio_name)) {
        colnames(Object@data)[colnames(Object@data) == Object@labeling_ratio_columns] = Object@labeling_ratio_name
      }
      else {
        colnames(Object@data)[colnames(Object@data) == Object@labeling_ratio_columns] = "LR"
        Object@labeling_ratio_name = "LR"
      }
    }
    if(!length(Object@labeling_ratio_columns)) {
      labeling_ratio_name = "LR"
      if(length(Object@labeling_ratio_name)) {
        labeling_ratio_name = Object@labeling_ratio_name
      } else {
        Object@labeling_ratio_name = labeling_ratio_name
      }
      if(trace) {
        cat(crayon::blue(paste0("Computing labeling ratio and inserting it into column `", labeling_ratio_name, "`.\n")))
      }
      LR = Object@data[[Object@intensity_name]]/(Object@data[[Object@intensity_name]] + Object@data[[Object@intensity_columns[1]]])
      if(Object@as_percentage) {
        LR = LR * 100
      }
      Object@data[,labeling_ratio_name] <- LR
    }
    
    if(!is.na(score_threshold)) {
      if(Object@higher_score_better) {
        Object@data = Object@data[(get(Object@score_name) >= score_threshold) & (get(Object@score_columns[1]) >= score_threshold)]
      } else {
        Object@data = Object@data[(get(Object@score_name) <= score_threshold) & (get(Object@score_columns[1]) <= score_threshold)]
      }
    }
    if(!is.na(labeling_ratio_threshold)) {
      Object@data = Object@data[get(Object@labeling_ratio_name) >= labeling_ratio_threshold]
    }
    Object@data$N <- 1
  }
  Object
}





