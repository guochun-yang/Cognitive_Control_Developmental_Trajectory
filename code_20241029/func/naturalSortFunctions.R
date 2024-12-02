# Custom sorting function with variable patterns
natural_sort <- function(files, prefix, suffix) {
  # Create a regular expression pattern based on the provided prefix and suffix
  pattern <- paste0("(?<=", prefix, ")[[:digit:]]+(?=", suffix, ")")
  
  # Extract the numeric parts of the filenames based on the pattern
  parts <- regmatches(files, gregexpr(pattern, files, perl = TRUE))
  
  # Convert the numeric parts to integers for proper comparison
  numeric_parts <- lapply(parts, function(x) as.integer(ifelse(length(x) > 0, x[[1]], NA)))
  
  # Order the files based on the numeric parts, ignoring files that don't match the pattern
  ordered_indices <- order(unlist(numeric_parts))
  
  # Return the sorted files
  return(files[ordered_indices])
}