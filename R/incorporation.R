cluster_rates <- function(Object,
                          var,
                          n,
                          timepoints = Object@timepoints,
                          cluster_method = "ward.D",
                          distance_method = "euclidean",
                          cluster_names = LETTERS[1:n],
                          sort = T)
{
  if(length(cluster_names) < n) {
    stop(paste("length of `cluster_names` must equal", n))
  }
  var_cols <- get_cols(Object, var, timepoints)
  mat <- impute(Object, vars = var, timepoints = timepoints)@master_tbl
  mat_var1 <- mat[,..var_cols]
  if(distance_method == "pearson" || distance_method == "spearman") {
    mat_dist <- as.dist(1 - cor(t(as.matrix(mat_var1)), method = distance_method))
  } else {
    mat_dist <- dist(mat_var1, method = distance_method)
  }
  hc <- hclust(mat_dist, method = cluster_method)
  k <- as.factor(cluster_names[cutree(hc, n)])
  mat_var2 <- mat_var1[, lapply(.SD, median), by = k]
  if(sort)
  {
    ord <- mat_var2[order(rowSums(mat_var2[, mapply(function(x, y) x - median(y), .SD[,-"k"], mat_var1)])), k]
    levels(k) <- setNames(lapply(ord, as.character), cluster_names[1:n])
  }
  k
}

heatmap <- function(Object,
                    var,
                    sample = NULL,
                    timepoints = Object@timepoints,
                    cluster_names = NULL,
                    split = NULL,
                    cluster_method = "ward.D",
                    distance_method = "euclidean",
                    filename = NULL,
                    width = NA,
                    height = NA,
                    sort = T,
                    plot = T,
                    column_title = Object@time_unit,
                    ...)
{
  if(is.null(split) & !is.null(cluster_names)) {
    split <- length(cluster_names)
  }
  cond = setdiff(names(Object@design), Object@time_unit)
  obj = Object
  if(length(sample) && length(cond)) {
    if(length(cond) == 1) {
      obj@master_tbl$test = obj@master_tbl[[cond]] == sample
    } else {
      obj@master_tbl = obj@master_tbl[, test := rowSums(.SD %in% sample, na.rm = T) == length(cond), .SDcols = cond]
    }
    for(i in seq_along(obj@consts)) {
      obj@consts[[i]] = obj@consts[[i]][obj@master_tbl$test]
    }
    obj@master_tbl = obj@master_tbl[(test),-"test"]
  }
  mat <- impute(obj, var = var, timepoints)@master_tbl[,get_cols(obj, var, timepoints),with=F]
  mat <- as.matrix(mat)
  colnames(mat) <- as.character(timepoints)
  q <- seq(0.25,0.75,length.out = 5)
  val <- sapply(q, quantile, x = mat)
  at <- c(seq(0, val[1], length.out = 4)[-4], val, seq(val[5], 100, length.out = 4)[-1])
  col_fun = circlize::colorRamp2(at, rev(RColorBrewer::brewer.pal(11, "Spectral")))
  dhm <- ComplexHeatmap::densityHeatmap(
    column_title = "",
    col = RColorBrewer::brewer.pal(9, "Blues"),
    mat, ylab = paste0(var, " (%)"),
    ylim = c(0,100),
    height = grid::unit(0.93, "npc"),
    heatmap_legend_param = list(
      title_position = "topcenter",
      grid_width = grid::unit(0.3, "in"),
      legend_direction = "horizontal",
      title = "Density"
    )
  )
  nam <- if(var == paste0("Heavy ", Object@incorporation_name)) {
    Object@incorporation_full_name
  } else {
    Object@labeling_ratio_full_name
  }
  hm <- ComplexHeatmap::Heatmap(
    mat,
    name = nam,
    col = col_fun,
    row_split = split,
    cluster_row_slices = T,
    cluster_rows = T,
    cluster_columns = F,
    show_row_names = F,
    row_title_rot = 0,
    clustering_distance_rows = distance_method,
    clustering_method_rows = cluster_method,
    na_col = "black",
    heatmap_legend_param = list(
      title_position = "leftcenter-rot",
      grid_width = grid::unit(2, "in"),
      col_fun = col_fun,
      legend_direction = "horizontal",
      side = "left"
    ),
    height = grid::unit(5, "cm"),
    column_names_rot = 0
  )
  ord <- ComplexHeatmap::row_order(hm)
  panel_fun = function(index, nm) {
    xrange = c(obj@time_zero, max(timepoints))
    xaxis = ComplexHeatmap::annotation_axis_grob(at = round(seq(xrange[1], xrange[2], length.out = 5)[-1]),
                                 side = "bottom", facing = "outside",
                                 gp = grid::gpar(cex = 0.6666667))
    yaxis = ComplexHeatmap::annotation_axis_grob(at = seq(0, 100, 25),
                                 side = "left", facing = "outside",
                                 gp = grid::gpar(cex = 0.6666667))
    grid::pushViewport(grid::viewport(xscale = xrange, yscale = c(0, 100)))
    grid::grid.rect(gp = grid::gpar(col = "black", fill = NA))
    grid::grid.draw(xaxis)
    grid::grid.draw(yaxis)
    cols <- get_cols(obj, var, timepoints)
    df <- obj@master_tbl[index, ..cols]
    x = seq(xrange[1], xrange[2], length.out = 100)
    cols <- paste0(rep(var, 100), " (", obj@time_unit, " ", x, ")")
    df2 <- data.table::as.data.table(model(obj, subset = index, var = var, x = x))
    colnames(df2) <- cols
    df[,paste0(var, " (", obj@time_unit, " ", 0, ")")] <- 0
    suppressWarnings(df <- data.table::melt(df))
    df <- df[, lapply(.SD, median, na.rm = T), by = "variable"]
    df$variable <- as.numeric(stringi::stri_extract_last_regex(df$variable, "\\d+"))
    df2[,paste0(var, " (", obj@time_unit, " ", 0, ")")] <- 0
    suppressWarnings(df2 <- data.table::melt(df2))
    df2 <- df2[, lapply(.SD, median, na.rm = T), by = "variable"]
    df2$variable <- as.numeric(stringi::stri_extract_last_regex(df2$variable, "[\\d\\.]+"))
    g1 <- ggplot2::ggplot(df, ggplot2::aes(x = variable, y = value)) +
      ggplot2::coord_cartesian(xlim = xrange, ylim = c(0,100)) +
      ggplot2::theme_void() + ggplot2::geom_point() + ggplot2::geom_line() +
      ggplot2::scale_x_continuous(expand = c(0,0), position = "top") +
      ggplot2::scale_y_continuous(expand = c(0,0), position = "right") +
      ggplot2::theme(legend.position = "none")
    
    grid::grid.draw(ggplot2::ggplotGrob(g1))
    grid::popViewport()
  }
  yaxis = function(index, nm) {
    grid::pushViewport(grid::viewport())
    g = grid::grid.text(label = "RITZ", x = grid::unit(0, "npc"),
                        vjust = 1,
                        gp = grid::gpar(fontsize = 8), rot = 90)
    grid::grid.draw(g)
    grid::popViewport()
  }
  anno = ComplexHeatmap::anno_zoom(
    align_to = ord, which = "row", panel_fun = panel_fun,
    link_width = grid::unit(6, "mm"),
    size = grid::unit(2, "cm"), width = grid::unit(4, "cm"), gap = grid::unit(5, "mm"),
    link_gp = grid::gpar(col = "grey", fill = NA), extend = grid::unit(0, "cm")
  )
  text_anno = ComplexHeatmap::anno_zoom(
    align_to = ord, which = "row", panel_fun = yaxis,
    link_width = grid::unit(6, "mm"),
    size = grid::unit(2, "cm"), width = grid::unit(0.4, "cm"),
    gap = grid::unit(0.5, "mm"), link_gp = grid::gpar(col = NA, fill = NA),
    extend = grid::unit(0, "cm")
  )
  lhs = rev(floor(lengths(ord)/2))
  rhs = rev(ceiling(lengths(ord)/2) - 1)
  words = NULL
  for(i in seq_len(4)) {
    words = c(words, c(rep("", lhs[i]), var, rep("", rhs[i])))
  }
  hm <- ComplexHeatmap::Heatmap(
    mat,
    right_annotation = ComplexHeatmap::rowAnnotation(
      foo = anno, 
      bar = text_anno
    ),
    name = nam,
    col = col_fun,
    row_split = split,
    cluster_row_slices = T,
    cluster_rows = T,
    cluster_columns = F,
    show_row_names = F,
    row_title_rot = 0,
    clustering_distance_rows = distance_method,
    clustering_method_rows = cluster_method,
    na_col = "black",
    heatmap_legend_param = list(
      legend_direction = "horizontal",
      title_position = "topcenter",
      grid_width = grid::unit(0.1, "in"),
      col_fun = col_fun
    ),
    height = grid::unit(1, "npc"),
    column_names_rot = 0,
    column_names_centered = T
  )
  if(!is.null(cluster_names)){
    rnk <- 1:split
    if(sort) {
      rnk <- rank(rowSums(sapply(as.data.frame(mat), function(x) sapply(ord, function(y) median(x[y]) - median(x)))))
    }
    hm@row_title = cluster_names[rnk]
  }
  lst <- ComplexHeatmap::`%v%`(dhm, hm)
  if(plot) {
    ComplexHeatmap::draw(lst, 
                         column_title = column_title,
                         column_title_side = "bottom",
                         ...)
  }
  if(!is.null(filename)) {
    save_plot(
      filename = filename, plot = lst, 
      width = width, height = height, 
      column_gap = grid::unit(0, "cm"),
      column_title = column_title,
      column_title_side = "bottom",
      ...
    )
  }
  lst
}

get_enrichment_from_cluster <- function(Object,
                                        group,
                                        var,
                                        n,
                                        sample = NULL,
                                        group_delimiter = ";",
                                        timepoints = Object@timepoints,
                                        cluster_names = LETTERS[1:n],
                                        cluster_method = "ward.D",
                                        distance_method = "euclidean",
                                        threshold = 0.05, 
                                        strategy = c("BH", "rFDR", "cFDR", "p-value"),
                                        sort = T) {
  
  strategy <- match.arg(strategy)
  cond = setdiff(names(Object@design), Object@time_unit)
  obj = Object
  if(length(sample) && length(cond)) {
    if(length(cond) == 1) {
      obj@master_tbl$test = obj@master_tbl[[cond]] == sample
    } else {
      obj@master_tbl = obj@master_tbl[, test := rowSums(.SD %in% sample, na.rm = T) == length(cond), .SDcols = cond]
    }
    for(i in seq_along(obj@consts)) {
      obj@consts[[i]] = obj@consts[[i]][obj@master_tbl$test]
    }
    obj@master_tbl = obj@master_tbl[(test),-"test"]
  }
  data <- impute(obj, vars = var, timepoints)@master_tbl
  data$Cluster <- cluster_rates(obj,
                                n = n,
                                timepoints = timepoints,
                                cluster_method = cluster_method,
                                distance_method = distance_method,
                                var = var,
                                cluster_names = cluster_names, sort = sort)
  f <- function(x) {
    unlist(strsplit(x, group_delimiter))
  }
  cols <- setdiff(names(data), group)
  data <- data[, lapply(.SD, f), by = cols, .SDcols = group]
  data[get(group) == "NA"][,group] = NA
  if(any(grepl("category", tolower(group))))
  {
    data <- data[!grepl("function unknown", tolower(get(group[grepl("category", tolower(group))])))]
    data <- data[!grepl("predict", tolower(get(group[grepl("category", tolower(group))])))]
    # data$`COG category` <- name2letter[data$`COG category`, letter, on = "name"]
  }
  data2 <- data[, .(group = data.table::melt(.SD, measure.vars = group)$variable,
                    name = data.table::melt(.SD, measure.vars = group)$value,
                    Cluster = Cluster), by = 1:nrow(data), .SDcols = group]
  mat <- data.table::dcast(data2, name + group ~ Cluster, fun = length, value.var = c("group"))
  mat <- na.omit(mat, c("group","name"))
  mat <- mat[rowSums(is.na(mat[,.(group, name)])) == 0,]
  cols <- colnames(mat)[-(1:2)]
  mat$r_total <- rowSums(mat[,..cols])
  mat$all_total <- nrow(data)
  data3 <- data[,.(N = .N), by = Cluster]
  mat <- cbind(mat, data.table::as.data.table(setNames(as.list(data3$N), paste0(data3$Cluster, "_tot"))))
  p_mat <- sapply(seq_along(cols), function(i) {
    all_total <- mat$all_total
    r_total <- mat$r_total
    c_total <- mat[[paste0(cols[i], "_tot")]]
    count <- mat[[cols[i]]]
    phyper(count - 1, r_total, all_total - r_total, c_total, lower.tail = F)
  })
  p <- sort(p_mat)
  if(strategy == "rFDR") {
    sig <- sapply(rank(p_mat, ties.method = "first"), function(i) {
      idx <- round(i * 1.6)
      if(idx > length(p)) {
        return(1)
      }
      fdr <- (p[idx]*length(p))/sum(p <= p[idx])
      if(fdr > 1) {
        return(1)
      }
      fdr
    })
  } else if(strategy == "cFDR") {
    sig <- sapply(rank(p_mat, ties.method = "first"), function(i) {
      c <- 0
      for(n in 1:i) {
        c <- c + 1/(i-n+1)
      }
      fdr <- c*((p[i]*length(p))/sum(p <= p[i]))
      if(fdr > 1) {
        return(1)
      }
      fdr
    })
  } else if (tolower(strategy) == "BH") {
    sig <- p.adjust(p_mat, method = "BH")
  } else {
    sig <- p_mat
  }
  sig_mat <- matrix(sig, ncol = ncol(p_mat))
  keep <- rowSums(sig_mat < threshold, na.rm = T) > 0
  if(all(!keep)) {
    warning(paste0("No enrichment in \"", group, "\"."))
    return(NULL)
  }
  sig_mat = sig_mat[keep,,drop=F]
  sig_mat[sig_mat >= threshold] <- NA_real_
  colnames(sig_mat) <- cols
  mat = mat[keep,1:2,with=F]
  sig_mat <- cbind(mat, sig_mat)
  sig_mat
}


plot_curves <- function(Object, vars, by = NULL, on = NULL, name = NULL, add_points = T,
                      top_n_by = NULL, min_n_by = NULL, top_n_on = NULL, min_n_on = NULL,
                      nrow = NULL, ncol = NULL, timepoints = Object@timepoints,
                      filename = NULL, plot = T, width = NA, height = NA, separate_vars = F,
                      tag_facets = F,tag_pool = LETTERS, tag_hjust = -0.5, tag_vjust = 1.5, tag_size = 3,
                      aggregate = c("median", "mean", "shrink", "ws", "weighted"), trace = T) {
  df <- Object@master_tbl
  grp <- unique(c(on, by))
  cols1 <- get_cols(Object, vars, timepoints)
  df <- na.omit(df, grp)
  if(is.null(name)){
    name <- unique(df[,..grp])
  }
  df <- df[name, , on = grp]
  if(length(on)) {
    keep_on <- df[, .N, by = on]
    keep_on <- keep_on[order(N)]
    top_n_on <- ifelse(is.null(top_n_on), ifelse(is.null(min_n_on), nrow(keep_on), sum(keep_on$N >= min_n_on)), top_n_on)
    keep_on <- keep_on[,tail(.SD, top_n_on)]
    df <- df[keep_on, on = on]
  }
  if(length(by) && grepl("category", by) ) {
    df <- df[tolower(get(by)) != "function unknown"]
  }
  if(length(on) && grepl("category", on) ) {
    df <- df[tolower(get(on)) != "function unknown"]
  }
  if(length(by)) {
    keep_by <- df[, .N, by = by]
    keep_by <- keep_by[order(N)]
    top_n_by <- ifelse(is.null(top_n_by), ifelse(is.null(min_n_by), nrow(keep_by), sum(keep_by$N >= min_n_by)), top_n_by)
    keep_by <- keep_by[,tail(.SD, top_n_by)]
    df <- df[keep_by, on = by]
  }
  if(length(by)) {
    df <- df[, setNames(list(unlist(strsplit(get(by), ";"))), by), by = setdiff(names(df), by)]
  }
  if(length(on)) {
    df <- df[, setNames(list(unlist(strsplit(get(on), ";"))), on), by = setdiff(names(df), on)]
  }
  aggregate = match.arg(aggregate)
  df <- aggregate_by(Object, data = df, vars = vars, by = grp, method = aggregate,
                     timepoints = timepoints, trace = trace)
  df2 <- df[,-..cols1]
  if(!is.list(timepoints)) {
    timepoints = list(timepoints)
  }
  if(length(timepoints) == 1 & length(vars) > 1) {
    timepoints = rep(timepoints, length(vars))
  }
  
  time_range = lapply(timepoints, function(t) range(c(Object@time_zero, t)))
  x = lapply(time_range, function(t) seq(t[1],t[2],length.out = 100))
  if(add_points) {
    x = mapply(function(y, z) sort(unique(c(y, z))), timepoints, x, SIMPLIFY = F)
  } else {
    mid_point = lapply(time_range, function(t) diff(t)/2)
    x = mapply(function(y, t) sort(unique(c(y, t))), mid_point, x, SIMPLIFY = F)
  }
  for(v in seq_along(vars)) {
    con = Object@consts[[vars[v]]]
    upper = unlist(unique(con[, grepl("_upper$", names(con)), with = F]))
    names(upper) = gsub("_upper$", "", names(upper))
    lower = unlist(unique(con[, grepl("_lower$", names(con)), with = F]))
    names(lower) = gsub("_lower$", "", names(lower))
    tmp = df[,get_cols(Object, vars[[v]], timepoints[[v]]),with=F]
    tmp = curve_fitting(Object, vars[v], timepoints = timepoints[[v]],
                        lower = lower, upper = upper, data = tmp)
    df2[, get_cols(Object, vars[[v]], x[[v]])] = data.table::as.data.table(
      model(x = x[[v]], var = vars[[v]], consts = tmp)
    )
  }
  df <- data.table::melt(df, id.vars = grp)
  df$x <- as.numeric(stringi::stri_extract_last_regex(df$variable, "\\d+"))
  df$variable <- stringi::stri_extract_first_regex(df$variable, paste0(vars, collapse = "|"))
  df2 <- data.table::melt(df2, id.vars = grp, value.name = "value1")
  df2$x <- as.numeric(stringi::stri_extract_last_regex(df2$variable, "[\\d\\.]+"))
  df2$variable <- stringi::stri_extract_first_regex(df2$variable, paste0(vars, collapse = "|"))
  cols3 <- setdiff(colnames(df2), "value1")
  if(!add_points) {
    df <- data.table::setnames(df2[,c(cols3, "value1"),with = F], c(cols3, "value"))
    ind <- unlist(mapply(function(m,v) {
      df[, .I[variable == v & x == mid_point]] 
    }, mid_point, vars, SIMPLIFY = F))
    df <- df[ind,]
  }
  df2 <- df[df2,, on = cols3]
  df2$variable2 = df2$variable
  if(length(on)) {
    df2$variable2 = df2[[on]]
  }
  df2$variable = paste0(df2$variable, " (%)")
  
  col = NULL
  n <- length(unique(df2$variable2))
  if(length(on)) {
    col = c('#222222', '#F3C300', '#875692', '#F38400', '#A1CAF1', '#BE0032', '#C2B280', '#848482', '#008856', '#E68FAC', '#0067A5', '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26')
    col = col[1:n]
    df2[["variable2"]] <- factor(df2[["variable2"]], levels = names(sort(table(df2[["variable2"]]), decreasing = F)))
  }
  g <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = value1, group = variable2)) +
    ggplot2::theme_classic() +
    # labs(title = stri_replace_first_fixed(name, " ", "\n")) +
    ggplot2::geom_line(eval(parse(text = ifelse(
      is.null(col),
      "ggplot2::aes(linetype = variable)",
      "ggplot2::aes(color = variable2)"
    )))) + 
    ggplot2::xlab("Time (day)") +
    ggplot2::geom_point(eval(parse(text = ifelse(
      is.null(col),
      "ggplot2::aes(x = x, y = value, shape = variable)",
      "ggplot2::aes(x = x, y = value, shape = variable2, color = variable2)"
    )))) +
    ggplot2::scale_shape_manual(values=1:n) +
    ggplot2::scale_color_manual(values = col) +
    ggplot2::coord_cartesian(ylim = c(0,100)) +
    ggplot2::guides(color = ggplot2::guide_legend(title = on),
           linetype = ggplot2::guide_legend(title = "y-axis"),
           shape = ggplot2::guide_legend(title = ifelse(!length(on), "y-axis", on)))
  if(length(by)) {
    by2 = by
    by2[!(by2 %in% make.names(by2))] <- paste0('`', by2[!(by2 %in% make.names(by2))], '`')
    labeller = parse(text = paste0("ggplot2::labeller(", by2, " = ggplot2::label_wrap_gen(30))"))
    g <- g + ggplot2::facet_wrap(eval(parse(text=paste0("ggplot2::vars(", by2,")"))),
                                 nrow = nrow, ncol = ncol, scales = "free", strip.position = "top",
                                 labeller = eval(labeller))
  } else if(separate_vars) {
    g <- g + ggplot2::facet_wrap(ggplot2::vars(variable), nrow = nrow, ncol = ncol,
                                 scales = "free", strip.position = "left")
  }
  g <- g + ggplot2::theme(strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(size = 12),
                     text = ggplot2::element_text(family = "sans"),
                     strip.placement = "outside",
                     axis.title.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.3, size = 12))
  if((length(by) || separate_vars) && tag_facets) {
    gb <- ggplot2::ggplot_build(g)
    lay <- gb$layout$layout
    tags <- cbind(lay, label = tag_pool[lay$PANEL], x = -Inf, y = Inf)
    g <- g + ggplot2::geom_text(data = tags, ggplot2::aes_string(x = "x", y = "y", label = "label"),
                                hjust = tag_hjust, inherit.aes = FALSE, vjust = tag_vjust, size = tag_size)
  }
   if(!is.null(filename)) {
    val = ifelse(is.null(col), max(strwidth(df2$variable, units = "inches", cex = circlize::fontsize(10))), max(strwidth(df2$variable2, units = "inches", cex = circlize::fontsize(10))))
    if(is.na(height) & is.na(width)) {
      n_panels <- length(unique(ggplot2::ggplot_build(g)$data[[1]]$PANEL))
      dims = ggplot2::wrap_dims(n_panels)
      width = 3 * dims[2]
      height = 2.25 * dims[1]
    }
    save_plot(filename = filename, plot = g, width = width + val, height = height)
  }
  if(plot) {
    grid::grid.draw(g)
  }
  g
}

plot_density <- function(Object, var, plot = T, filename = NULL,
                              width = NA, height = NA,
                              ylim = c(0,100), ylab = paste0(var, " (%)"),
                              title = "", column_names_rot = 0,
                              legend_title_position = "topcenter",
                              legend_direction = "horizontal",
                              legend_title = "Density",
                              heatmap_legend_side = "top",
                              ...
) {
  mat <- setnames(Object@master_tbl[,get_cols(Object, var),with=F], as.character(Object@timepoints))
  den <- ComplexHeatmap::densityHeatmap(
    mat, ylim = ylim, ylab = ylab, title = title,
    column_names_rot = column_names_rot,
    heatmap_legend_param = list(
      title_position = legend_title_position,
      legend_direction = legend_direction,
      title = legend_title
    )
  )
  if(plot) {
    ComplexHeatmap::draw(den, heatmap_legend_side = heatmap_legend_side, ...)
  }
  if(length(filename)) {
    save_plot(filename = filename, plot = den, width = width, height = height, column_gap = grid::unit(0, "cm"), 
              heatmap_legend_side = heatmap_legend_side, ...)
  }
  den
}
