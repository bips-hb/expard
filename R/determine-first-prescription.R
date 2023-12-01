#' First Prescription
#' 
#' \code{determine_first_perscription} returns a time point
#' for the first prescription of the drug of interest 
#' given the time frame of the patient (\code{simulation_time})
#' and the chance for it be prescribed is \code{min_chance_drug}.
#' The patient, in this case, is always prescribed the drug at least
#' once. The time point follows a \eqn{\Gamma} distribution. 
#' 
#' \emph{Note:} we assume here that the probability of 
#' the drug being prescribed does not depend on any previous prescriptions. 
#' 
#' @param n_patients Number of patients (Default: 1)
#' @param simulation_time The total number of time points for this patient
#' @param min_chance_drug The probability of the drug being prescribed 
#'                        at any given time point
#'                        
#' @return A time point between \code{1} and \code{simulation_time}
#' @examples 
#' set.seed(1)
#' determine_first_prescription(n_patients = 4, simulation_time = 10, min_chance_drug = 0.1)
#' # -> [1] 2 3 5 9 
#' @export
determine_first_prescription <- function(n_patients, simulation_time, min_chance_drug) { 
  
  # check correctness input --------
  if (n_patients < 1)      { stop("n_patients should be >= 1") }
  if (simulation_time < 1) { stop("simulation_time should be >= 1") }
  if (min_chance_drug <= 0 || min_chance_drug > 1) { 
    stop("min_chance_drug should be in the interval (0,1]") 
  }
  
  # determining the probabilities that the drug is prescribed at each
  # time point (1, 2, ..., simulation_time). Follows a Gamma distribution
  probs <- sapply(1:simulation_time, function(t) (1 - min_chance_drug)^(t - 1))
  probs <- probs * min_chance_drug / sum(probs) # normalization
  
  # sample a random time point given the probabilities
  return(
    sample(x = 1:simulation_time, n_patients, replace = TRUE, prob = probs)
  )
}