# Minimal JSON writer in base R (the pipeline's DESeq2 container has no jsonlite). Sourced by host_verify_de.R and
# host_verify_gsea.R. Lists with names become objects, unnamed lists arrays, data frames arrays of row objects;
# NA, NaN and Inf become null.
json_esc <- function(s) gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", s))

json_scalar <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (length(v) == 0 || is.na(v)) return("null")
  if (is.logical(v)) return(if (v) "true" else "false")
  if (is.numeric(v)) {
    if (!is.finite(v)) return("null")
    if (v == round(v) && abs(v) < 1e15) return(sprintf("%.0f", v))
    return(formatC(v, digits = 6, format = "g"))
  }
  paste0('"', json_esc(as.character(v)), '"')
}

to_json <- function(x, ind = 0) {
  p0 <- strrep("  ", ind); p1 <- strrep("  ", ind + 1)
  if (is.null(x)) return("null")
  if (is.data.frame(x)) {
    x <- lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
    if (length(x) == 0) return("[]")
  }
  if (is.list(x)) {
    if (length(x) == 0) return("[]")
    v <- vapply(x, to_json, "", ind = ind + 1)
    nm <- names(x)
    if (!is.null(nm) && all(nzchar(nm)))
      return(paste0("{\n", paste0(p1, '"', json_esc(nm), '": ', v, collapse = ",\n"), "\n", p0, "}"))
    return(paste0("[\n", paste0(p1, v, collapse = ",\n"), "\n", p0, "]"))
  }
  if (is.factor(x)) x <- as.character(x)
  if (length(x) == 1 && is.null(names(x))) return(json_scalar(x))
  v <- vapply(seq_along(x), function(i) json_scalar(x[[i]]), "")
  nm <- names(x)
  if (!is.null(nm) && all(nzchar(nm)))
    return(paste0("{", paste0('"', json_esc(nm), '": ', v, collapse = ", "), "}"))
  paste0("[", paste(v, collapse = ", "), "]")
}
