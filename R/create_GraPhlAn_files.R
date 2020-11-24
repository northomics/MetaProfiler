create_GraPhlAn_files <- function(Object, graphlan_dir, var, day, filename, output_dir = getwd(), ranks = get_ranks(Object@master_tbl), by = ranks[1], min_distinct_peptide = 2, min_cell_nb = 2, node_size = 50) {
  value_tbl <- NULL
  by_is_rank <- any(ranks %in% by)
  if (by_is_rank & ranks[1] != by) {
    stop("`by` seem to be a phylogenetic rank. This is only be valid if ranks start at `by` (e.g. rank = c(\"Superkingdom\", \"Kingdom\", \"Phylum\", \"Class\", \"Order\", \"Family\", \"Genus\", \"Species\"), by = \"Superkingdom\").")
  }
  Object@annotation_table[grepl("Lacto", Object@annotation_table$species)]
  data <- na.omit(data.table::copy(Object@master_tbl), unique(by, ranks[1]))

  if(grepl("category", by)) {
    data = data[get(by) != "Function unknown"]
    data = data[, setNames(list(unlist(strsplit(get(by), ";"))),by), by = setdiff(colnames(data), by)]
  }
  lca_column = colnames(data)[grepl("lca", tolower(colnames(data)))]
  rank_column = colnames(data)[grepl("rank", tolower(colnames(data)))]
  for(r in ranks) {
    by2 <- c(r, by)
    cols <- c(lca_column, by, rank_column, "value", "p.value")
    data2 <- data.table::copy(data)
    data2 <- data2[, test := (.N >= min_distinct_peptide) & (length(na.omit(unlist(.SD))) >= 2), by = by2, .SDcols = get_cols(Object, var, day)]
    data2 <- data2[(test), -"test"]
    data2 <- na.omit(data2[,.(r = r, t.test(na.omit(unlist(.SD)))$estimate, t.test(na.omit(unlist(.SD)))$p.value), by = by2,  .SDcols = get_cols(Object, var, day)])
    value_tbl <- rbind(value_tbl, data.table::setnames(data2, cols), fill = T)
  }
  value_tbl$p.value = p.adjust(value_tbl$p.value, "fdr")
  value_tbl = value_tbl[p.value < 0.05]
  check_cols = c(lca_column, rank_column, by)
  keep = value_tbl[unique(data[,..check_cols]),, on = check_cols]
  keep = keep[!is.na(keep$value)]
  cols <- unique(c(by, ranks))
  all_lineage <- unique(data[unique(keep[,..check_cols]), ..cols, on = check_cols])
  unique_lineage <- all_lineage
  all_lineage <- unique(data.table::rbindlist(apply(all_lineage, 1, function(x) {
    x <- zoo::na.trim(x)
    x[is.na(x)] <- ""
    res <- NULL
    for(i in 1:length(x)) {
      tmp <- x[1:i]
      if(tmp[i] == "") {
        tmp[tmp == ""] <- NA 
      }
      tax <- as.list(setNames(c(tmp, rep(NA, length(x)-i)), names(x)))
      res <- c(res, list(tax))
    }
    data.table::rbindlist(res, fill = T)
  }), fill = T))
  # all_lineage <- cbind(organism = "root", all_lineage)
  
  unique_lineage <- unique(unique_lineage[order(rowSums(!is.na(unique_lineage)), decreasing = T),])
  phylo_tmp <- list(unique_lineage[1,])
  for(i in 2:nrow(unique_lineage))
  {
    is_unique <- T
    for(j in 1:length(phylo_tmp))
    {
      if(all(na.omit(unlist(unique_lineage[i,])) %in% na.omit(unlist(phylo_tmp[[j]]))))
      {
        is_unique <- F
      }
    }
    if(is_unique) {
      phylo_tmp <- c(phylo_tmp, list(as.list(unique_lineage[i,])))
    }
  }
  
  unique_lineage <- data.table::rbindlist(phylo_tmp)
  
  
  
  # input <- apply(unique_lineage, 1, function(x) {
  #   x <- x[!is.na(x)]
  #   embed(rev(x), 2)
  # })
  # el <- unique(as.data.table(do.call(rbind, input)))
  # el <- el[, lapply(.SD, function(x) gsub("\\W", "_", x))]
  # levels <- unique(na.omit(unlist(all_lineage)))
  # require(igraph)
  # g <- graph_from_data_frame(as.data.frame(el), directed = F)
  # nwk <- graph_to_newick(g, "root")
  # tree <- read.tree(text = nwk)
  # d <- data.frame(node= (1:(Ntip(tree) + Nnode(tree)))[-(Ntip(tree) + 1)], images=paste0(gsub("^\\w__", "", c(tree$tip.label, tree$node.label[-1])), ".png"))
  # pg <- ggtree(tree) %<+% d + xlim(NA, 9)
  # pg <- pg + geom_nodelab(aes(image=images), geom="image", size = .09) +
  #   geom_tiplab(aes(image=images), geom="image", size = .09, offset=2, align=T, hjust=0)
  # grid.draw(pg)
  # save_plot("mygraph4.pdf", pg, width = 8.27, height = 11.69)
  # require(ape)
  # identify(tree)
  # 
  # tree$node.label[21 - Ntip(tree)]
  # Ntip(tree)
  # names(V(g))
  # require(png)
  # imgs <- sapply(paste0(gsub("^\\w__", "", c(tree$tip.label, tree$node.label)), ".png"), readPNG)
  # V(g)$raster <- imgs
  # require(ggraph)
  # layout <- layout.reingold.tilford(g)
  # layout <- -layout[,2:1]   
  # png("mygraph.png", height=30, width=60, units = "in", res = 300)
  # plot(g, vertex.shape = "raster", vertex.label=NA,
  #           layout = layout,
  #           vertex.size = 10, vertex.size2 = 10*5/3,
  #           asp = 0.35, margin = -0.1)
  # dev.off()
  # list(edge = as.matrix(el), Nnodes = Nnodes, tip.label = tip.label)
  # plot(g, vertex.shape = "raster", vertex.label=NA, layout = layout_nicely(g))
  # attr(tree, "class") <- "all_lineage"
  
  input <- paste0(apply(all_lineage, 1, function(x) {
    x <- x[!is.na(x)]
    paste(paste0(tolower(substring(names(x), 1, 1)), "__", gsub("\\W+", "_", x), collapse = "."))
  }), collapse = "\n")
  n <- length(na.omit(unique(all_lineage[[by]])))
  col = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
  annot_back <- data.table::setnames(data.table::data.table(col[1:n], na.omit(unique(all_lineage[[by]]))), c("col", by))
  q <- seq(0.25,0.75,length.out = 5)
  val <- sapply(q, quantile, x = unique(value_tbl$value))
  at <- c(seq(0, val[1], length.out = 4)[-4], val, seq(val[5], 100, length.out = 4)[-1])
  col_fun = circlize::colorRamp2(at, rev(RColorBrewer::brewer.pal(11, "Spectral")))
  header <- paste("title_font_size\t13",
                  "start_rotation\t90",
                  "annotation_background_alpha\t0.15",
                  "clade_separation\t0.35",
                  "annotation_legend_font_size\t5",
                  "annotation_font_size\t6",
                  "annotation_font_stretch\t0",
                  "clade_marker_size\t3",
                  "branch_bracket_depth\t0.5",
                  "branch_thickness\t1.2", sep = "\n")
  annot <- paste0(apply(all_lineage, 1, function(x) {
    x <- x[!is.na(x)]
    if(!by_is_rank & length(x) == 1) {
      return(paste0(
        paste(paste0(tolower(substring(names(x), 1, 1)), "__", gsub("\\W+", "_", x), collapse = "."),
              "clade_marker_size\t0", sep = "\t"),
        "\n",
        paste(paste0(tolower(substring(names(x), 1, 1)), "__", gsub("\\W+", "_", x), collapse = "."),
              "clade_marker_color\t#000000", sep = "\t")
      ))
    }
    on = c(rank_column,lca_column, by)
    to_get <- data.table::setnames(data.table::data.table(x[length(x)], x[by], names(x)[length(x)]), c(lca_column, by, rank_column))
    val <- value_tbl[to_get, value, on = on]
    paste0(
      paste(paste0(tolower(substring(names(x), 1, 1)), "__", gsub("\\W+", "_", x), collapse = "."),
            "clade_marker_size", node_size, sep = "\t"),
      "\n",
      paste(paste0(tolower(substring(names(x), 1, 1)), "__", gsub("\\W+", "_", x), collapse = "."),
            "clade_marker_color", col_fun(val), sep = "\t")
    )
  }), collapse = "\n")
  id <- NULL
  tax <- NULL
  b <- NULL
  rnk <- NULL
  for(r in ranks) {
    phylo_tmp <- data.table::copy(unique_lineage)
    by2 <- unique(c(r, by))
    phylo_tmp <- phylo_tmp[, test := .N < min_cell_nb, by = by2]
    nm <- na.omit(unique(phylo_tmp[(test), ..by2]))
    nm$tmp <- ifelse(grepl(" ", nm[[r]]), paste0(tolower(substring(r,1,1)),
                                               substring(nm[[r]],1,1),
                                               stringi::stri_trim_both(stringi::stri_extract_first_regex(nm[[r]], " \\w"))), 
                     paste0(tolower(substring(r,1,1)), substring(nm[[r]],1,1)))
    convert_to_unique <- unique(nm[, c(r, "tmp"), with = F])
    convert_to_unique$tmp <- make.unique(convert_to_unique$tmp, sep = "")
    nm$tmp <- convert_to_unique[nm[[r]], tmp, on = r]
    tax <- c(tax, nm[[r]])
    id <- c(id, nm$tmp)
    rnk <- c(rnk, rep(r, length(nm[[r]])))
    b <- c(b, nm[[by]])
  }
  tbl2 <- unique(na.omit(data.table::setnames(data.table::data.table(rnk, tax, b, id), c(rank_column, lca_column, by, "l"))))
  if(grepl("category", by)) {
    Letter_NOG_Category <- data.table::fread("Letter_NOG_Category.csv")
    cat <- unique(value_tbl[[by]])
    tbl2 <- rbind(data.table::setnames(data.table::data.table(by, cat, cat, Letter_NOG_Category[cat, Letter, on = "Name"]), c(rank_column, lca_column, by, "l")), tbl2)
  }
  
  annot <- paste(annot, paste0(apply(all_lineage, 1, function(x) {
    x <- x[!is.na(x)]
    t <- tbl2[data.table::setnames(data.table::data.table(x[length(x)], x[1], names(x)[length(x)]), c(lca_column, by, rank_column)), on = c(lca_column, by, rank_column)]
    if(!is.na(t$l)) {
      paste(paste0(tolower(substring(names(x), 1, 1)), "__", gsub("\\W+", "_", x), collapse = "."),
            "\tannotation\t", t$l, sep = "")
    } else {
      paste(paste0(tolower(substring(names(x), 1, 1)), "__", gsub("\\W+", "_", x), collapse = "."),
            "annotation", t[[lca_column]], sep = "\t")
    }
  }), collapse = "\n"), sep = "\n")
  annot <- paste(annot, paste0(apply(all_lineage, 1, function(x) {
    x <- x[!is.na(x)]
    if(!by_is_rank & (length(x) == 1)) {
      return(paste(paste0(tolower(substring(names(x), 1, 1)), "__", gsub("\\W+", "_", x), collapse = "."),
                   "annotation_background_color\t#FFFFFF", sep = "\t"))
    } 
    col <- annot_back[x[1], col, on = by]
    paste(paste0(tolower(substring(names(x), 1, 1)), "__", gsub("\\W+", "_", x), collapse = "."),
          "annotation_background_color", col, sep = "\t")
  }), collapse = "\n"), sep = "\n")
  
  annot <- paste(header, annot, sep = "\n")         
  
  write(annot, file = "annot.txt")
  
  write(input, file = "input.txt")
  path = gsub("(\\w):/", "\\L\\1/", output_dir, perl = T)
  command = c(paste("export PATH=", file.path(graphlan_dir , ":$PATH"), sep=''),
              paste0("graphlan_annotate.py \"/mnt/", path, "/input.txt\" \"/mnt/",
                     path, "/plot.xml\" --annot \"/mnt/", path, "/annot.txt\""),
              paste0("graphlan.py \"/mnt/", path, "/plot.xml\" \"/mnt/",
                     path, "/plot.png\" --dpi 300 --size 5 --pad 0"))
  writeLines("Copy and paste these lines in your ubuntu subsystem:\n")
  writeLines(command)
  readline("\nOnce the files are generated, press enter to continue...")
  
  dat <- unique(tbl2[, c("l", lca_column, rank_column), with = F])
  dat <- dat[, .SD[order(l)], by = rank_column]
  at <- apply(dat[dat[[rank_column]] %in% ranks, c("l", lca_column), with = F], 1, paste0, collapse = ": ")
  
  annot_lgnd <- ComplexHeatmap::Legend(at, title = "Annotation Legend:", legend_gp = grid::gpar(fill = NA, fontsize = 12, fontfamily = "sans"), ncol = 4)
  annot_width = as.numeric(grid::convertUnit(grid::widthDetails(annot_lgnd@grob), "in")) + 0.08
  annot_height = as.numeric(grid::convertUnit(grid::heightDetails(annot_lgnd@grob), "in"))
  annot_height = annot_height + annot_height/ceiling(length(at)/4) + 0.08
  
  save_plot("annotation_legend.png", annot_lgnd, width = annot_width, height = annot_height)
  color_key = ComplexHeatmap::Legend(col_fun = col_fun, title = paste0(Object@labeling_ratio_full_name, " (%)"), title_position = "leftcenter-rot", grid_height = grid::unit(10, "mm"))
  color_key_width = as.numeric(grid::convertUnit(grid::widthDetails(color_key@grob), "in")) + 0.08
  color_key_height = as.numeric(grid::convertUnit(grid::heightDetails(color_key@grob), "in")) + 0.08
  save_plot(filename = "LR_color_key.png", color_key, width = color_key_width, height = color_key_height)
  mp <- grid::rasterGrob(png::readPNG("plot.png"))
  al <- grid::rasterGrob(png::readPNG("annotation_legend.png"))
  ck <- grid::rasterGrob(png::readPNG("LR_color_key.png"))
  
  if(grepl("category", by)) {
    cog_lgnd <- paste0(Letter_NOG_Category[annot_back[[by]], Letter, on = "Name"], ": ", annot_back[[by]])
    cog_lgnd <- ComplexHeatmap::Legend(cog_lgnd[order(nchar(cog_lgnd))], title = "COG Category Legend:", legend_gp = grid::gpar(fill = annot_back$col[order(nchar(cog_lgnd))]), ncol = 3)
    cog_width = as.numeric(grid::convertUnit(grid::widthDetails(cog_lgnd@grob), "in")) + 0.1
    cog_height = as.numeric(grid::convertUnit(grid::heightDetails(cog_lgnd@grob), "in"))
    cog_height = cog_height + cog_height/ceiling(length(at)/3) + 0.08
    save_plot("col2nog.png", cog_lgnd, width = cog_width, height = cog_height)
    cl <- grid::rasterGrob(png::readPNG("col2nog.png"))
    ag <- gridExtra::arrangeGrob(grobs = list(cl, mp, al, ck), layout_matrix = matrix(c(1,2,2,2,3,NA,4,NA,NA,NA), ncol = 2), heights = c(cog_height * annot_width/cog_width, rep(annot_width/3, 3), annot_height), widths = c(annot_width, color_key_width * (annot_width/3)/color_key_height))
    save_plot(filename, ag, width = sum(c(annot_width, color_key_width * (annot_width/3)/color_key_height)), height = sum(c(cog_height * annot_width/cog_width, annot_width, annot_height)))
  } else {
    phyla_lgnd <- annot_back[[by]]
    phyla_lgnd <- ComplexHeatmap::Legend(phyla_lgnd[order(nchar(phyla_lgnd))], title = "Phyla Legend:", legend_gp = grid::gpar(fill = annot_back$col[order(nchar(phyla_lgnd))]), ncol = 4)
    phyla_width = as.numeric(grid::convertUnit(grid::widthDetails(phyla_lgnd@grob), "in")) + 0.1
    phyla_height = as.numeric(grid::convertUnit(grid::heightDetails(phyla_lgnd@grob), "in"))
    phyla_height = phyla_height + phyla_height/ceiling(length(at)/3) + 0.08
    save_plot("col2phyla.png", phyla_lgnd, width = phyla_width, height = phyla_height)
    cl <- grid::rasterGrob(png::readPNG("col2phyla.png"))
    ag <- gridExtra::arrangeGrob(grobs = list(cl, mp, al, ck), layout_matrix = matrix(c(1,2,2,2,3,NA,4,NA,NA,NA), ncol = 2), heights = c(phyla_height * annot_width/phyla_width, rep(annot_width/3, 3), annot_height), widths = c(annot_width, color_key_width * (annot_width/3)/color_key_height))
    save_plot(filename, ag, width = sum(c(annot_width, color_key_width * (annot_width/3)/color_key_height)), height = sum(c(annot_width, annot_height)))
  }
}

