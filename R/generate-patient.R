#' Generate a Patient 
#' 
#' Generates a patient, which is a collection of different 
#' drug-ADR pairs, see function \code{\link{generate_drug_ADR_pair}}. 
#' 
#' @param simulation_time The total number of time steps
#' @param n_drug_ADR_pairs Number of drug-ADR pairs simulated
#' @param risk_model Vector with risk models. Each risk model is given as 
#'                   a string, e.g., \code{"risk_model_current_use()"}. The 
#'                   vector must have the length \code{n_drug_ADR_pairs}#' 
#' @param min_chance_drug Vector with the probabilities of the drug being prescribed
#'            when the drug history and the ADR history have no effect. 
#'            Must have a length of \code{n_drug_ADR_pairs}. 
#' @param avg_duration Average number of time points a patient is exposed
#'                     once exposed to the drug (Determines \code{max_chance})
#'                     Must have a length of \code{n_drug_ADR_pairs}. 
#' @param max_chance_drug The probability of the ADR when 
#'            the drug history and the ADR history have the highest 
#'            possible effect (Default: \code{NULL})
#'            Must have a length of \code{n_drug_ADR_pairs}. 
#' @param guaranteed_exposed If \code{TRUE}, the patient is exposed to the drug
#'            at least once
#'            Must have a length of \code{n_drug_ADR_pairs}. 
#' @param min_chance The probability of the ADR when 
#'            the drug history has no effect
#'            Must have a length of \code{n_drug_ADR_pairs}. 
#' @param max_chance The probability of the ADR when 
#'            the drug history has the highest possible effect
#'            Must have a length of \code{n_drug_ADR_pairs}. 
#' 
#' @return A list of drug-ADR pairs, see \code{\link{generate_drug_ADR_pair}}
#'
#' @seealso \code{\link{generate_cohort}}, \code{\link{generate_drug_ADR_pair}}, \code{\link{generate_drug_history}}  
#' @export
generate_patient <- function(simulation_time = 100,
                             n_drug_ADR_pairs = 50,
                             risk_model = rep("risk_model_current_use()", n_drug_ADR_pairs),
                             min_chance_drug = rep(.1, n_drug_ADR_pairs),
                             avg_duration = rep(5, n_drug_ADR_pairs),
                             max_chance_drug = rep(NULL, n_drug_ADR_pairs),
                             guaranteed_exposed = rep(TRUE, n_drug_ADR_pairs),
                             min_chance  = rep(.1, n_drug_ADR_pairs),
                             max_chance  = rep(.4, n_drug_ADR_pairs)) { 
  
  res <- lapply(1:n_drug_ADR_pairs, function(i) { 
    
    generate_drug_ADR_pair(simulation_time = simulation_time, 
                           eval( parse(text = risk_model[i]) ),
                           min_chance_drug[i], 
                           avg_duration[i], 
                           max_chance_drug[i],
                           guaranteed_exposed[i], 
                           min_chance[i], 
                           max_chance[i])
    })
  
  class(res) <- c(class(res), "patient")
  return(res)
}

#' Function for printing a patient 
#' @export
print.patient <- function(patient) { 
  cat(sprintf("Patient\n")) 
  
  simulation_time 
  
  cat(sprintf("\nNo. of time points: %d\n\n", patient$simulation_time))
  
  lapply(1:length(patient), function(i) { 
    cat(sprintf("Drug-ADR pair %d\n", i))
    print(patient[[i]]) 
  })
  
  cat(sprintf("\n"))
}
