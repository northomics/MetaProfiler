unipept_api = function(input, type = "pept2funct", base = "http://api.unipept.ugent.be/api/v1/", 
                       extra = T, equate_il = T, domains = F) {
  
  input_len = length(input)
  entries = NULL
  while(length(input)) {
    call = paste0(base, type, ".json?input[]=", paste0(input, collapse = "&input[]="), ifelse(equate_il, "&equate_il=true", ""),
                  ifelse(extra, "&extra=true", ""), ifelse(domains, "&domains=true", ""))
    i = NULL
    while(nchar(call) >= 2000) {
      if(!length(i)) {
        i = 1
      } else {
        i = c(i, i[length(i)] + 1)
      }
      call = paste0(base, type, ".json?input[]=", paste0(input[-i], collapse = "&input[]="),
                    ifelse(equate_il, "&equate_il=true", ""),
                    ifelse(extra, "&extra=true", ""), ifelse(domains, "&domains=true", ""))
    }
    json = httr::GET(call)
    entry = data.table::as.data.table(jsonlite::fromJSON(content(json, "text"), flatten = T))
    ids = as.character(seq_len(nrow(entry)))
    fun = function(x) {
      if(!is.list(x)) return(x)
      x = data.table::rbindlist(x, fill = T, idcol = "id")
      x = x[, lapply(.SD[max(protein_count, na.rm = T) == protein_count], paste0, collapse = ";"), by = id]
      x$id = as.character(x$id)
      x = x[ids, , on = "id"][,-"id"]
      names(x)[1] = gsub("\\w+_", "", names(x)[1])
      x
    }
    entry = entry[, do.call(cbind, lapply(.SD, fun))]
    names(entry) = gsub("\\.", "_", names(entry))
    entries = rbind(entries, entry)
    if(!length(i)) input = NULL
    else input = input[i]
    cat(paste0(input_len - length(input), " out of ", input_len))
  }
  entries
}
