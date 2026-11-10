# Writes a data frame as CSV. Each double is written with the fewest significant digits
# (15 to 17) that read back as the same number, so the processed tables round-trip exactly.
write_table <- function(x, path) {
  text <- vapply(x, is.character, logical(1))
  x[] <- lapply(x, function(v) if (is.double(v)) shortest_digits(v) else v)
  utils::write.csv(x, path, row.names = FALSE, na = "NA", quote = which(text))
  invisible(path)
}

shortest_digits <- function(v) {
  out <- formatC(v, digits = 15, format = "g")
  known <- !is.na(v)
  for (digits in 16:17) {
    off <- known
    off[known] <- as.numeric(out[known]) != v[known]
    if (!any(off)) break
    out[off] <- formatC(v[off], digits = digits, format = "g")
  }
  out <- trimws(out)
  out[is.na(v)] <- NA
  out
}

read_table <- function(path) utils::read.csv(path, stringsAsFactors = FALSE)
