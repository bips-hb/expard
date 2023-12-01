#' Generate Drug Prescription History
#' 
#' Returns whether or not the drug is prescribed 
#' at the next time point (\code{1}) or not (\code{0}). 
#' The probability depends on the drug prescriptions, 
#' the ADR history and the \code{drug_model}.
#' 
#' @param simulation_time Total number of time points
#' @param min_chance The probability of the drug being prescribed
#'            when the drug history and the ADR history have no effect
#' @param avg_duration Average number of time points a patient is exposed
#'                     once exposed to the drug (Determines \code{max_chance})
#' @param max_chance The probability of the ADR when 
#'            the drug history and the ADR history have the highest 
#'            possible effect (Default: \code{NULL})
#' @param guaranteed_exposed If \code{TRUE}, the patient is exposed to the drug
#'            at least once
#' 
#' @return \code{1} or \code{0}
#' @family Update functions
#' @examples 
#' generate_drug_history(simulation_time = 10, 
#'                       min_chance = .1,
#'                       avg_duration = 5
#'                       guaranteed_exposed = TRUE)
#' @export
generate_drug_history <- function(simulation_time = 10, 
                                  min_chance,
                                  avg_duration = 5, 
                                  max_chance = NULL, 
                                  guaranteed_exposed = TRUE,
                                  ...) { 
  
  if (avg_duration < 1) { 
    stop("avg duration must be at least 1")  
  }
  
  if (is.null(max_chance)) { 
    max_chance = (avg_duration - 1) / avg_duration 
  }
  
  # determine first time the drug is prescribed by a random sample
  if (guaranteed_exposed) { 
    first <- determine_first_prescription(n_patients = 1, simulation_time, min_chance) 
    
    if (first == simulation_time) { 
      drug_history <- c(rep(0, simulation_time - 1), 1) 
      return(drug_history)
    }
    
    # initialize the drug and ADR time processes 
    drug_history <- c(rep(0, first - 1), 1, rep(NA, simulation_time - first))
    
    # simulate the remaining time points
    for (t in (first+1):simulation_time) { 
      # get the probability of the next time point
      prob <- min_chance + (max_chance - min_chance) * drug_history[t-1]
      drug_history[t] <- rbinom(1,1,prob)
    }
    
    
  } else { 
    
    # determine the first time point
    drug_history <- c(rbinom(1,1, min_chance), rep(NA, simulation_time - 1))
    
    # simulate the remaining time points
    for (t in 2:simulation_time) { 
      # get the probability of the next time point
      prob <- min_chance + (max_chance - min_chance) * drug_history[t-1]
      drug_history[t] <- rbinom(1,1,prob)
    }
  }
  
  return(drug_history)
}
