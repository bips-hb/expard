#' Compute Frequency Table
#' 
#' Determine how often different values from time_steps_ago 
#' appear in the data set and how often the ADR occurs for 
#' each of them. This is used for speeding up the computation of the 
#' loglikeihood
#' 
#' @param mat Matrix with relevant values (depends on model)
#' @param adr_history History of ADRs
#' 
#' @return Tibble
#' @export
determine_frequency_unique_values <- function(mat, adr_history) {
  unique_values <- sort(unique(as.vector(mat)))
  
  freq_table <-
    lapply(sort(unique_values), function(unique_value) {
      # where does this 'unique_value' occur in the data
      indices <- which(mat == unique_value)
      
      # count how often an ADR occurred at those time points
      n_adr <- sum(adr_history[indices])
      
      c(
        unique_value = unique_value,
        freq       = length(indices),
        # how often does it occur
        n_adr      = n_adr,
        n_no_adr   = length(indices) - n_adr
      )
    })
  
  # create a tibble and return it
  return(dplyr::as_tibble(do.call(rbind, freq_table)))
}