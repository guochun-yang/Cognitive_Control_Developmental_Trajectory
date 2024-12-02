# Function to perform FDR correction
fdr_correction <- function(pmtx_adjust) {
  return(p.adjust(p, method = p.adjust.methods, n = length(p)))
}