#' \eqn{2 \times 2} Tables
#' 
#' Creates the \eqn{2 \times 2} contingency tables for *all* 
#' drug-ADR pairs in a specific cohort. The function is basically
#' a wrapper for \code{\link{create2x2table}}
#' 
#' \cr\cr
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
#' @param cohort A cohort; see \code{\link{generate_cohort}} 
#' @param method Method used to construct the table; either 
#'               \code{time-point}, \code{drug-era} and \code{patient}. 
#'               See the description for more information (Default: 
#'               \code{time-point})
#'
#' @return A list of \code{cont_table} objects
#'
#' @references 
#' Zorych, I., Madigan, D., Ryan, P., & Bate, A. (2013). Disproportionality methods for 
#' pharmacovigilance in longitudinal observational databases. 
#' Statistical Methods in Medical Research, 22(1), 39–56. 
#' https://doi.org/10.1177/0962280211403602  
#' @seealso \code{\link{create2x2table}}
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
create2x2tables <- function(cohort,
                            method = c("time-point",
                                       "drug-era",
                                       "patient"),
                            verbose = TRUE) {
  
  
  if (!(method[1] %in% c("time-point", "drug-era", "patient"))) { 
    stop(sprintf("method should be either '%s', '%s' or '%s'", 
                 "time-point", "drug-era", "patient")) 
  }
  
  if (verbose) { 
    cat("Generating 2x2 tables...\n")
    pb <- txtProgressBar(min = 0, max = cohort$n_drug_ADR_pairs, style = 3) 
  }
  
  # initialize tables
  tables <- lapply(1:cohort$n_drug_ADR_pairs, function(i) { 
    table <- create2x2table(cohort[[i]], method)
    
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
    
    return(table)
  })
  
  if (verbose) {
    close(pb) 
    cat("DONE generating tables...\n")
  }
  
  return(tables)
}
