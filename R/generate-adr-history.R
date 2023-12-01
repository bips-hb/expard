#' Returns ADR History
#' 
#' Generates a binary vector of ADR occurrences given a drug history
#' and a risk model. 
#' 
#' @param drug_history Binary vector denoting the 
#'            drug prescriptions history
#' @param risk_model A risk model, e.g., \code{\link{risk_model_current_use}}
#' @param min_chance The probability of the ADR when 
#'            the drug history has no effect 
#' @param max_chance The probability of the ADR when 
#'            the drug history has the highest possible effect
#' 
#' @return Binary vector of the same length as \code{drug_history}
#' @export
generate_adr_history <- function(drug_history, 
                                 risk_model, 
                                 min_chance,
                                 max_chance, 
                                 ...) { 
  
  simulation_time <- length(drug_history)

  prob <- min_chance + (max_chance - min_chance) * risk_model(drug_history, min_chance, max_chance)
  
  adr_history <- sapply(prob, function(p) rbinom(1,1,p))
  return(adr_history)
}