#' Generate a Drug-ADR pair
#' 
#' Generates the drug- and ADR history
#' 
#' @param simulation_time The total number of time steps
#' @param risk_model One of the risk models
#' @param min_chance_drug The probability of the drug being prescribed
#'            when the drug history and the ADR history have no effect
#' @param avg_duration Average number of time points a patient is exposed
#'                     once exposed to the drug (Determines \code{max_chance})
#' @param max_chance_drug The probability of the ADR when 
#'            the drug history and the ADR history have the highest 
#'            possible effect (Default: \code{NULL})
#' @param guaranteed_exposed If \code{TRUE}, the patient is exposed to the drug
#'            at least once
#' @param min_chance The probability of the ADR when 
#'            the drug history has no effect 
#' @param max_chance The probability of the ADR when 
#'            the drug history has the highest possible effect
#' 
#' @return A list with
#'       \item{\code{drug_history}}{binary vector of the 
#'                drug history}
#'       \item{\code{adr_history}}{binary vector of the 
#'                ADR history}
#'
#' @seealso \code{\link{generate_cohort}} 
#' @examples 
#' drug_ADR_pair <- expard::generate_drug_ADR_pair()
#' @export
generate_drug_ADR_pair <- function(simulation_time = 100,
                                   risk_model = risk_model_current_use(),
                                   min_chance_drug = .1,
                                   avg_duration = 5,
                                   max_chance_drug = NULL,
                                   guaranteed_exposed = TRUE,
                                   min_chance  = .1,
                                   max_chance  = .4)
{ 
  
  
  drug_history <- generate_drug_history(simulation_time, min_chance_drug, 
                                        avg_duration, max_chance_drug, 
                                        guaranteed_exposed)
  
  
  adr_history <- generate_adr_history(drug_history, risk_model, min_chance, max_chance) 
  
  res <- list(
    drug_history = drug_history,
    adr_history = adr_history
  )
  
  class(res) <- c(class(res), "drug_ADR_pair")
  
  return(res)
}

#' Function for printing a drug-ADR pair
#' @export
print.drug_ADR_pair <- function(drug_ADR_pair) { 

  simulation_time <- length(drug_ADR_pair$drug_history)
  
  cat(sprintf("Drug: "))
  sapply(drug_ADR_pair$drug_history, function(x) { 
    if (x == 1) { 
      cat(green(1))  
    } else {
      cat(blue("."))
    }
  })
  
  cat(sprintf("\nADR:  "))
  sapply(drug_ADR_pair$adr_history, function(y) { 
    if (y == 1) { 
      cat(red(1))  
    } else {
      cat(blue("."))
    }
  })
  
  cat(sprintf("\n"))
}
