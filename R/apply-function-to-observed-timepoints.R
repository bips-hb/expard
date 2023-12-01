#' Apply Function to Subset Drug History
#' 
#' Observations might not be present for certain patients, since they entered 
#' the data set later time = 1 or before the end. This function applies, for 
#' example, a risk model only to the observed time points. 
#' 
#' @param drug_history The drug history. Unobserved time points are denoted by 
#'                     \code{NA}
#' @param fn A function (Default \code{\link{risk_model_current_use}}). 
#' 
#' @return Vector of the same length as \code{drug_history}. None observed time points
#'         are denoted by \code{NA}
#' @export                 
apply_function_to_observed_timepoints <- function(drug_history, 
                                                  fn = expard::risk_model_current_use()) { 
  
  simulation_time <- length(drug_history)
  
  # determine the indices that are not NA
  indices_not_NA <- which(!is.na(drug_history))
  
  # get only the observed drug history given these indices
  observed_drug_history <- drug_history[indices_not_NA]
  
  # determine the values for only this part of the drug history
  values_fn <- fn(observed_drug_history)
  
  # initialize the risk vector with the same length as the original drug history
  values <- rep(NA, simulation_time)
  
  # fill in the observed risks
  values[indices_not_NA] <- values_fn
  
  return(values)
}