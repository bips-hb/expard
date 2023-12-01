#' Creating a \eqn{2 \times 2} Table
#' 
#' Creates a \eqn{2 \times 2} contingency table for a 
#' drug and ADR. \cr\cr
#' A table is structured in the following form:  
#' \tabular{lcc}{
#'    \tab ADR \tab not ADR\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#' 
#' # Ways to construct the table 
#' The counts can be constructed in three different ways. See for 
#' two of them Zorych et al. (2013).  
#' 
#' ## By individual time-points (\code{time-point}) 
#' Each time-point is counted separately. The counts are the 
#' \emph{number of time points} that
#' \itemize{ 
#'    \item{\code{a} - the drug was prescribed and the ADR occurred}
#'    \item{\code{b} - the drug was not prescribed but the ADR occurred}
#'    \item{\code{c} - the drug was prescribed but the ADR did not occur}
#'    \item{\code{d} - the drug was not prescribed and the ADR did not occur}
#' }
#' Note that in this case the total count \code{n = a + b + c + d} is 
#' the same as the total number of time points observed, i.e., the total 
#' number of patients times the number of time points observed: 
#' \code{n_patients * simulation_time}. 
#' 
#' ## By individual patients (\code{patient})
#' In this case, individual patients are counted. The counts are the \emph{number 
#' of patients} that 
#' \itemize{
#'    \item{\code{a} - were prescribed the drug and did experience the ADR }
#'    \item{\code{b} - were never prescribed the drug and did experience the ADR}
#'    \item{\code{c} - were prescribed the drug and never experienced the ADR}
#'    \item{\code{d} - were never prescribed the drug and never experienced the ADR}
#' }
#' In this case, the total count \code{n = a + b + c +d} is the same as 
#' the number of patients, \code{n_patients}.
#' 
#' ## By individual patients (\code{drug-era})
#' In this case we look at \emph{drug-eras}, i.e., periods in which the 
#' patients was prescribed or not prescribed the drug for a longer time. 
#' For example, if the patient was prescribed the drug from time point 3 to 
#' 6, then that period is called a drug-era. 
#' The counts are the \emph{number of drug- and non-drug eras} in which 
#' \itemize{ 
#'    \item{\code{a} - the drug was prescribed and the ADR occurred}
#'    \item{\code{b} - the drug was not prescribed but the ADR occurred}
#'    \item{\code{c} - the drug was prescribed but the ADR did not occur}
#'    \item{\code{d} - the drug was not prescribed and the ADR did not occur}
#' }
#' In this case, the total count \code{n} is the total number of drug- and 
#' non-drug eras.
#' 
#' @param drug_ADR_pair A drug-ADR pair; see \code{\link{generate_cohort}} 
#' @param method Method used to construct the table; either 
#'               \code{time-point}, \code{drug-era} and \code{patient}. 
#'               See the description for more information (Default: 
#'               \code{time-point})
#'
#' @return A \code{cont_table} object; a list with  
#'     \item{\code{a},\code{b},\code{c},\code{d}}{The counts in the table}
#'     \item{\code{method}}{Method that was used to construct the table}
#'     \item{\code{n}}{Total count (\code{n = a + b + c + d}). The interpretation 
#'                depends on the \code{method} used} 
#'
#' @references 
#' Zorych, I., Madigan, D., Ryan, P., & Bate, A. (2013). Disproportionality methods for 
#' pharmacovigilance in longitudinal observational databases. 
#' Statistical Methods in Medical Research, 22(1), 39–56. 
#' https://doi.org/10.1177/0962280211403602  
#' @seealso \code{\link{generate_cohort}}
#' @examples 
#' set.seed(1)
#' cohort <- generate_cohort(n_patients = 200) 
#' 
#' # create the 2x2 contingency table per time-point, 
#' # drug-era and patient: 
#' create2x2table(cohort, method = "time-point")
#' create2x2table(cohort, method = "drug-era")
#' create2x2table(cohort, method = "patient")
#' @export
create2x2table <- function(drug_ADR_pair, method = c("time-point", 
                                                     "drug-era", 
                                                     "patient")) { 
  
  if (!(method[1] %in% c("time-point", "drug-era", "patient"))) { 
    stop(sprintf("method should be either '%s', '%s' or '%s'", 
                 "time-point", "drug-era", "patient")) 
  }
  
  # initialize tables
  table <- list(a = 0, b = 0, c = 0, d = 0, method = method[1])
  class(table) <- "cont_table"
    
  if (method[1] == "time-point") {
    
    drug <- drug_ADR_pair$drug_history == 1 
    ADR <- drug_ADR_pair$adr_history == 1
    
    table$a <- sum(drug & ADR, na.rm = TRUE)
    table$b <- sum(!drug & ADR, na.rm = TRUE)
    table$c <- sum(drug & !ADR, na.rm = TRUE)
    table$d <- sum(!drug & !ADR, na.rm = TRUE)
  }
  
  if (method[1] == "patient") {
    # go over all patients
    drug <- rowSums(as.matrix(drug_ADR_pair$drug_history), na.rm = TRUE) > 0
    ADR <- rowSums(as.matrix(drug_ADR_pair$adr_history), na.rm = TRUE) > 0
    
    table$a <- sum(drug & ADR, na.rm = TRUE)
    table$b <- sum(!drug & ADR, na.rm = TRUE)
    table$c <- sum(drug & !ADR, na.rm = TRUE)
    table$d <- sum(!drug & !ADR, na.rm = TRUE)
  }
  
  if (method[1] == "drug-era") { 
    
    
    
    sapply(1:cohort$n_patients, function(k) {
      
      # remove any none observed time points. They are represented by NAs
      indices_observed_time_points_drug <- which(!is.na(drug_ADR_pair$drug_history[k, ]))
      indices_observed_time_points_adr <- which(!is.na(drug_ADR_pair$adr_history[k, ]))
      indices_observed_time_points <- union(indices_observed_time_points_drug, indices_observed_time_points_adr)
      
      drug_history_patient <- drug_ADR_pair$drug_history[k, indices_observed_time_points]
      adr_history_patient <- drug_ADR_pair$adr_history[k, indices_observed_time_points]
    
      simulation_time_patient <- length(drug_history_patient)
      
      # first initialize some variables to keep track
      # in which era we (drug or non-drug) and whether
      # the ADR occured during this era
      in_drug_era <-
        drug_history_patient[1] == 1  # are we currently in a drug era?
      ADR_happened <-
        drug_history_patient[1] == 1    # did the ADR occur during this era?
      
      sapply(2:(simulation_time_patient - 1), function(t) {
        if (in_drug_era) {
          # in drug-era
          if (drug_history_patient[t]) {
            # drug prescribed on time point t?
            if (drug_history_patient[t]) {
              # did the ADR occur?
              ADR_happened <- TRUE
            }
          } else {
            # switch from a drug-era to a non-drug era
            in_drug_era <- FALSE
            if (ADR_happened) {
              table$a <<- table$a + 1
            } else {
              table$c <<- table$c + 1
            }
          }
        } else {
          # not in drug-era
          if (drug_history_patient[t] == 0) {
            # drug not prescribed
            if (drug_history_patient[t] == 1) {
              # ADR occurred
              ADR_happened <- TRUE
            }
          } else {
            # switch from a non-drug-era to a drug era
            in_drug_era <- TRUE
            if (ADR_happened) {
              table$b <<- table$b + 1
            } else {
              table$d <<- table$d + 1
            }
          }
          
        }
        ADR_happened <- adr_history_patient[t] == 1
      })
    })
  }
  
  
  table$n <- table$a + table$b + table$c + table$d 
  
  return(table)
}

#' Print function for 2x2 tables
#' @export
print.cont_table <- function(table) { 
  cat(sprintf("2 x 2 Contingency Table\n"))
  cat(sprintf("\tusing method '%s'\n\n", table$method))
  cat("         |\tADR\t|     not ADR\t| total\n")
  cat("------------------------------------------------\n")
  cat(sprintf("    drug |\t%d\t|\t%d\t| %d\n", table$a, table$c, table$a + table$c)) 
  cat(sprintf("not drug |\t%d\t|\t%d\t| %d\n", table$b, table$d, table$b + table$d)) 
  cat("------------------------------------------------\n")
  cat(sprintf("   total |\t%d\t|\t%d\t| %d\n", table$a + table$b, 
              table$c + table$d, 
              table$a + table$b + table$c + table$d)) 
  
  if ((table$a + table$c) == (table$a + table$b + table$c + table$d)) { 
    cat(crayon::magenta(sprintf("\nwarning: since the number of patients that were prescribed \nthe drug and the total number of patients is the same,\nit might be that the cohort was created like this on purpose"))) 
  }
}