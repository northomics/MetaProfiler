cluster_features <- function(Object, peptide_centric = T, peptide_charge_column = "guess",
                             id_type = c("both", "id", "feature"),
                             cluster_by = c(Object@incorporation_name, Object@labeling_ratio_name),
                             radius = c(0.1,0.1), distance_method = c("relative", "relative"),
                             merge_overlaps = T, element_wise = T, trace = T, progress = T) {
  if(!length(Object@peptide_column_no_PTMs) && !length(Object@peptide_column_PTMs) && !length(Object@accession_column)) {
    stop("At least one peptide column or protein column must be specified.")
  }
  if(peptide_centric) {
    if(length(peptide_charge_column) && peptide_charge_column == "guess") {
      charge = colnames(Object@data)[grepl("charge", tolower(colnames(Object@data)))]
      if(trace) {
        cat(crayon::silver("Using `", charge, "` as the peptide charge column.\n"))
      }
    }
    group_by = c(Object@peptide_column_no_PTMs, Object@peptide_column_PTMs, Object@accession_column, colnames(Object@design), charge)
  } else {
    group_by = c(Object@accession_column, colnames(Object@design))
  }
  if(progress) {
    writeLines("clustering features...")
  }
  centers = names(Object@data)[sapply(Object@data, is.numeric)]
  data = Object@data
  if(length(data$Feature)) {
    id_type = match.arg(id_type)
    if(id_type != "both") {
      data = data[Feature == id_type]
    }
  }
  clust <- qtclust(
    data,
    cluster_by = cluster_by,
    radius = radius,
    distance_method = distance_method,
    centers = centers,
    group_by = group_by,
    element_wise = element_wise,
    merge_overlaps = merge_overlaps,
    progress = progress
  )
  if(trace) {
    message(paste0("clustered ", nrow(Object@data), " features to ", nrow(clust$centers), "."))
  }
  clust$centers
  Object@data = clust$centers
  Object
}

qtclust <- function(data, cluster_by, radius, distance_method, group_by = NULL, centers = NULL, 
                    start = NULL, end = NULL, merge_overlaps = T, element_wise = F, progress = F) {
  if(is.null(centers))
  {
    centers = cluster_by
  }
  cols_needed <- unique(c(cluster_by, centers, group_by))
  gr <- data.table::copy(data[,..cols_needed])
  gr$idx <- 1:nrow(data)
  if(!is.null(group_by)) {
    gr <- gr[, id := .GRP, by = group_by]
    gr <- gr[order(id)]
  } else {
    gr$id <- 1
  }
  
  if(element_wise)
  {
    if(length(distance_method) != length(cluster_by) & length(distance_method) != 1) {
      stop(paste("length of distance_method must be either 1 or", length(cluster_by), "for clustering using element-wise distances"))
    }
    if(length(distance_method) == 1)
    {
      distance_method = rep(distance_method, length(cluster_by))
    }
    if (is.null(start) & is.null(end))
    {
      start = 1
      end = length(cluster_by)
    }
  } else {
    if (length(start) != length(distance_method) & length(end) != length(distance_method)) {
      stop("length of start and end must be the same as distance_method for clustering using row-wise distances")
    }
    if (is.null(start) & is.null(end))
    {
      start = 1
      end = length(cluster_by)
    }
  }
  require(Rcpp)
  sourceCpp("~MetaProfiler/qtclust.cpp")
  clust <- qtclust_c (as.matrix(gr[,..cluster_by]), n_groups = max(gr$id), id = gr$id,
                      groups = as.matrix(gr[,..centers]), radius = radius,
                      method = distance_method, start = start - 1, end = end - 1,
                      element_wise = element_wise, merge_overlaps = merge_overlaps,
                      verbose = progress)
  clust$centers <- data.table::setnames(data.table::as.data.table(clust$centers), c(centers, "N"))
  gr$cluster <- clust$cluster 
  clust$centers[,group_by] <- unique(gr[,c(group_by, "cluster"), with = F])[,-"cluster"]
  clust$SD <- data.table::setnames(data.table::as.data.table(clust$SD), centers)
  clust
}




