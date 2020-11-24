.accession_razor <- function(data, protein, peptide, accession_delimiter, progress) {
  x <- data.table::copy(data)
  if(is.list(x[[protein]])) {
    x$unique = as.numeric(lengths(x[[protein]]) == 1)
  } else {
    x <- x[, unique := as.numeric(!grepl(accession_delimiter, get(protein))) , by = protein]
    x[[protein]] <- stringi::stri_split_regex(x[[protein]], accession_delimiter)
  }
  x <- x[, lapply(.SD, unlist), by=1:nrow(x)][,-1]
  x <- x[, c("count", "unique") := .(.N, sum(unique)), by = protein]
  x <- x[order(count, unique, decreasing = T)]
  x <- x[, .(get(protein)[[1]], paste0(count, collapse = ","), paste0(unique, collapse = ",")), by = peptide]
  x <- data.table::setnames(x, c(peptide, paste(protein, "(razor)"), "count", "unique"))
  x$razor_count = as.numeric(stringi::stri_extract_first_regex(x$count, "\\d+"))
  x$razor_unique = as.numeric(stringi::stri_extract_first_regex(x$unique, "\\d+"))
  data <- x[data, , on = peptide]
  data
  # groups <- data.table::as.data.table(razor(x[[peptide]], x[[protein]], x[["unique"]]))
  # groups <- groups[, c(lapply(.SD, function(.sd) sapply(.sd, "[", 1)), list(x = x)), by = 1:nrow(groups), .SDcols = c("id", "count", "unique")]
  # groups <- groups[, .(x = unlist(x)), by = .(id, count, unique)]
  # groups <- groups[, order()]
  # groups[,peptide] <- x[[peptide]]
  # groups$Proteins
}

.make_annotation_table <- function(Object,
                                   compute_razor_protein,
                                   pep2pro,
                                   pep2func,
                                   pep2taxon,
                                   pro2func,
                                   trace,
                                   progress) {
  
  if(!any(colnames(Object@data) == Object@time_unit)) {
    warning(paste0("No column named `", Object@time_unit, "` in slot `data`. Please add/rename the column for the timepoints or change slot `time_unit` so that it matches the timepoint column in slot `data`."))
  }
  annotation_table <- data.table::data.table()
  if(length(pep2pro)) {
    delimiter <- table(unlist(stringi::stri_extract_all_regex(pep2pro[[Object@pep2pro_accession_column]], "[\t;,]")))
    delimiter <- names(delimiter)[which.max(delimiter)]
    if(!length(delimiter)) delimiter = ";"
    pep2pro = data.table::setnames(pep2pro[,paste0(unique(get(Object@pep2pro_accession_column)), collapse = delimiter), by = c(Object@pep2pro_peptide_column)], c(Object@pep2pro_peptide_column,Object@pep2pro_accession_column))
    if(compute_razor_protein) {
      if(trace) {
        cat(crayon::blue("Computing razor protein."))
      }
      pep2pro = .accession_razor(pep2pro, Object@pep2pro_accession_column, Object@pep2pro_peptide_column, delimiter, progress)      
    }
    peptides = unique(pep2pro[Object@data[[Object@annotate_with]], , on = Object@pep2pro_peptide_column] )
    annotation_table = data.table::data.table(Proteins = peptides[[Object@pep2pro_accession_column]], Peptides = peptides[[Object@pep2pro_peptide_column]])
    # annotation_table <- na.omit(annotation_table, "Proteins")
  } else if(length(Object@accession_column)) {
    peptides = unique(Object@data[, na.omit(c(Object@accession_column, Object@annotate_with)), with = F])
    annotation_table = data.table::data.table(Proteins = peptides[[Object@accession_column]], Peptides = peptides[[Object@annotate_with]])
    if(compute_razor_protein) {
      if(progress) {
        cat("Computing razor protein...")
      }
      annotation_table$Proteins = .accession_razor(annotation_table, "Proteins", "Peptides", "[\t;,]", progress)      
      if(progress) {
        cat("done.\n")
      }
    }
    annotation_table <- annotation_table[, .(Proteins = paste0(Proteins, collapse = ";")), by = "Peptides"]
    annotation_table <- unique(na.omit(annotation_table, "Proteins"))
  }
  if(!nrow(annotation_table)) {
    warning("It looks like the accession column is empty. Please provide the Peptides file that contains the assigned proteins or specify the accession column in the result file.")
  }
  if(length(Object@accession_column) && nrow(annotation_table)) {
    Object@data[,Object@accession_column] <- annotation_table[Object@data[[Object@annotate_with]], Proteins, on = "Peptides"]
  } else if (length(annotation_table$Proteins)) {
    Object@data[,"Proteins"] <- annotation_table[Object@data[[Object@annotate_with]], Proteins, on = "Peptides"]
    Object@accession_column = "Proteins"
  }
  if(length(pep2taxon)) {
    pep2taxon_columns = colnames(pep2taxon)
    lca_column = pep2taxon_columns[grepl("lca", tolower(pep2taxon_columns))]
    lca_rank_column = pep2taxon_columns[grepl("rank", tolower(pep2taxon_columns))]
    tax <- pep2taxon[annotation_table$Peptides, c(Object@pep2taxon_peptide_column, lca_column, lca_rank_column, Object@rank_columns), with = F, on = c(Object@pep2taxon_peptide_column)]
    names(tax)[1] = "Peptides"
    annotation_table <- tax[annotation_table, , on = "Peptides"]
  }
  if(length(pro2func) || length(pep2func)) {
    pro2func[pro2func == ""] <- NA
    if(!length(Object@pro2func_function_columns) && !length(Object@pep2func_function_columns)) {
      stop("Names for the functional annotation columns are empty. Please specify the columns.")
    }
    if(length(Object@pep2func_function_columns)) {
      function_columns = Object@pep2func_function_columns
      by = Object@pep2func_peptide_column
      search_by = "Peptides"
      delimiter <- table(unlist(stringi::stri_extract_all_regex(pep2func[[function_columns[[1]]]], "[\t;,]")))
      delimiter <- names(delimiter)[which.max(delimiter)]
      if(!length(delimiter)) delimiter = ";"
      func = data.table::setnames(pep2func[,lapply(.SD, function(x) paste0(unique(x), collapse = delimiter)), by = by, .SDcols = function_columns], c(by,function_columns))
    } else if(length(Object@pro2func_function_columns)) {
      function_columns = Object@pro2func_function_columns
      by = Object@pro2func_accession_column
      search_by = "Proteins"
      delimiter <- table(unlist(stringi::stri_extract_all_regex(pro2func[[function_columns[[1]]]], "[\t;,]")))
      delimiter <- names(delimiter)[which.max(delimiter)]
      if(!length(delimiter)) delimiter = ";"
      func = data.table::setnames(pro2func[,lapply(.SD, function(x) paste0(unique(x), collapse = delimiter)), by = by, .SDcols = function_columns], c(by,function_columns))
    }
    names(func)[names(func) %in% by] = search_by
    func = unique(func)
    annotation_table <- func[annotation_table,, on = c(search_by)]
  }
  Object@annotation_table = annotation_table
  Object
}

get_annotations <- function(Object, sample) {
  annotation_table <-  Object@annotation_table
  names(annotation_table)[names(annotation_table) == "Peptides"] = Object@annotate_with
  names(annotation_table)[names(annotation_table) == "Proteins"] = Object@accession_column
  on = c(Object@annotate_with, Object@accession_column, sample)
  Object@master_tbl <- annotation_table[Object@master_tbl, , on = c(on)]
  Object
}



