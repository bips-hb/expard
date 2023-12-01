#' Generate a Cohort 
#' 
#' Generates a cohort which is a collection of drug-ADR pair. 
#' It is a wrapper function for the \code{\link{generate_patient}} 
#' function. 
#' 
#' @param n_patients Number of patients 
#' @inheritParams generate_patient
#' @param verbose Show progress bar (Default: \code{FALSE})
#' 
#' @return A \code{cohort} object; A list of length \code{n_drug_ADR_pairs}. 
#'         Each entry contains two binary matrices; one for the drug histories
#'         of all patients and one for the ADR histories of all patients. Each 
#'         row is a patient
#' 
#' @seealso \code{\link{generate_patient}}   

#' @export
generate_cohort <- function(n_patients = 100, 
                            simulation_time = 100, 
                            n_drug_ADR_pairs = 50, 
                            risk_model = rep("risk_model_current_use()", n_drug_ADR_pairs), 
                            min_chance_drug = rep(.1, n_drug_ADR_pairs), 
                            avg_duration = rep(5, n_drug_ADR_pairs), 
                            max_chance_drug = rep(NULL, n_drug_ADR_pairs),
                            prob_guaranteed_exposed = rep(1, n_drug_ADR_pairs), 
                            min_chance  = rep(.1, n_drug_ADR_pairs), 
                            max_chance  = rep(.4, n_drug_ADR_pairs),
                            verbose = FALSE) { 
  
  
  if (verbose) { 
    cat("Generating the patients...\n")
    pb <- txtProgressBar(min = 0, max = n_patients, style = 3)
  }
  
  # generate patients 
  patients <- lapply(1:n_patients, function(i) { 
    
    guaranteed_exp <- rbinom(n_drug_ADR_pairs, 1, prob_guaranteed_exposed) 
    
    patient <- generate_patient(
      simulation_time,
      n_drug_ADR_pairs,
      risk_model,
      min_chance_drug,
      avg_duration,
      max_chance_drug,
      guaranteed_exposed = guaranteed_exp,
      min_chance,
      max_chance
    )
    
    if (verbose) { 
      setTxtProgressBar(pb, i)
    }
    return(patient)
  })
  
  if (verbose) { 
    close(pb) 
    cat("DONE generating the patients...\n\nOrganizing the data into matrices...\n")
    pb <- txtProgressBar(min = 0, max = n_drug_ADR_pairs, style = 3)
  }
  
  res <- lapply(1:n_drug_ADR_pairs, function(i) {
    drug_history <- matrix(rep(NA, simulation_time * n_patients), nrow = n_patients) 
    adr_history <- drug_history
    
    sapply(1:n_patients, function(p) {
      #print(patients[[p]][[i]])
      drug_history[p, ] <<- patients[[p]][[i]]$drug_history
      adr_history[p, ] <<- patients[[p]][[i]]$adr_history
    })

    if (verbose) { 
      setTxtProgressBar(pb, i)
    }
    
    list(drug_history = Matrix(drug_history, sparse = TRUE), 
         adr_history = Matrix(adr_history, sparse = TRUE))
  })

  if (verbose) { 
    close(pb) 
    cat("DONE organizing the data into matrices...\n")
  }
  
  res$n_patients <- n_patients
  res$n_drug_ADR_pairs <- n_drug_ADR_pairs
  res$simulation_time <- simulation_time
  
  class(res) <- "cohort"
  return(res)
}

#' Print function for \code{\link{generate_cohort}}
#' @export
print.cohort <- function(cohort) { 
  
  cat("Cohort\n\n")
  cat(sprintf("  No. patients:\t\t%d\n", cohort$n_patients))
  cat(sprintf("  No. drug-ADR-pairs:\t%d\n", cohort$n_drug_ADR_pairs))
  cat(sprintf("  No. time points:\t%d\n", cohort$simulation_time))
}