#' Negative of the Log-Likelihood Function
#' 
#' \code{loglikelihood} returns the negative 
#' log-likelihood function, defined as 
#' \deqn{\pi = \frac{\exp(\beta_0 + \code{risks}\cdot \beta)}{1 + \exp(\beta_0 + \code{risks}\cdot \beta)}} 
#' where \eqn{\pi} is the probability. 
#' 
#' @param betas Vector with two entries: \code{beta0, beta}
#' @param risks Risks determined by the risk model. The risk 
#'              model was applied to the observed drug prescription 
#'              history
#' @export
loglikelihood <- function(betas,
                          risks,
                          adr_history) {

  # extract parameters
  beta0 <- betas[1]
  beta <- betas[2]

  P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
  -1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
}


#' @export
loglikelihood_past <- function(param,
                               past,
                               drug_history,
                               adr_history) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  
  risk_model <- expard::risk_model_past(past)
  
  n_patients <- nrow(drug_history)
  simulation_time <- ncol(drug_history)
  
  # 'convert' the drug prescriptions. They reflect which period is considered
  # to have an increased risk
  risks <- matrix(0, nrow = n_patients, ncol = simulation_time)
  
  # go over all patients 
  for (i in 1:n_patients) { 
    # go over all timepoints 
    risks[i, ] <- risk_model(drug_history[i, ])
  }
  
  #risks <- Matrix(apply(drug_history, 1, function(d) risk_model(d)))
  
  P <- exp(beta0 + risks * beta) / (1 + exp(beta0 + risks * beta))
  -1 * sum( adr_history * log(P) + (1 - adr_history)*log(1 - P))
}


#' @export
loglikelihood_withdrawal <- function(param,
                            freq_table) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  phi <- param[3]
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  phi <- param[3]
  
  # transform to regular parameters
  pi0 <- exp(beta0) / (1 + exp(beta0))
  pi1 <- exp(beta) / (1 + exp(beta))
  rate <- exp(phi)
  
  # determine the risks given the time_since last exposure and the rate 
  # of the withdrawal model
  freq_table <- freq_table %>% mutate(
    risk_value = exp(-rate * (unique_value - 1))
  )
  
  # the risk when the time_since is 0 (never exposed or currently exposed)
  # is zero
  freq_table$risk_value[1] <- 0
  
  # determine the log-likelihood given the risks
  freq_table <- freq_table %>% rowwise() %>% 
    mutate(loglikelihood = -1*(n_adr*log((pi1 - pi0)*risk_value + pi0) + 
                                 n_no_adr * log(1 - (pi1 - pi0)*risk_value - pi0)))
  
  sum(freq_table$loglikelihood)
}



#' @export
loglikelihood_delayed <- function(param,
                                  freq_table) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  mu_log <- param[3]
  sigma_log <- param[4]
  
  # transform to regular parameters
  pi0 <- exp(beta0) / (1 + exp(beta0))
  pi1 <- exp(beta) / (1 + exp(beta))
  mu <- exp(mu_log)
  sigma <- exp(sigma_log)
  
  # to make sure that the highest value is indeed 1
  normalizing_factor <- dnorm(mu, mu, sigma)
  
  # determine the risks given the time_since last exposure and the rate 
  # of the withdrawal model
  freq_table <- freq_table %>% mutate(
    risk_value = dnorm(unique_value, mu, sigma) / normalizing_factor
  )
  
  # the risk when the time_since is 0 (never exposed or currently exposed)
  # is zero
  freq_table$risk_value[1] <- 0
  
  # determine the log-likelihood given the risks
  freq_table <- freq_table %>% rowwise() %>% 
    mutate(loglikelihood = -1*(n_adr*log((pi1 - pi0)*risk_value + pi0) + 
                                 n_no_adr * log(1 - (pi1 - pi0)*risk_value - pi0)))
  
  sum(freq_table$loglikelihood)
}





#' @export
loglikelihood_decaying <- function(param,
                                   freq_table) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  rate_log <- param[3]
  
  # transform to regular parameters
  pi0 <- exp(beta0) / (1 + exp(beta0))
  pi1 <- exp(beta) / (1 + exp(beta))
  rate <- exp(rate_log)
  
  # determine the risks given the time_since last exposure and the rate 
  # of the withdrawal model
  freq_table <- freq_table %>% mutate(
    #risk_value = exp(-rate * (unique_value - 1))
    risk_value = exp(-rate * (unique_value - 1))
  )
  
  # the risk when the time_since is 0 (never exposed or currently exposed)
  # is zero
  freq_table$risk_value[1] <- 0
  
  # determine the log-likelihood given the risks
  freq_table <- freq_table %>% rowwise() %>% 
    mutate(loglikelihood = -1*(n_adr*log((pi1 - pi0)*risk_value + pi0) + 
                                 n_no_adr * log(1 - (pi1 - pi0)*risk_value - pi0)))
  
  sum(freq_table$loglikelihood)
}





#' @export
loglikelihood_delayed_decaying <- function(param,
                                  freq_table) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  mu_log <- param[3]
  sigma_log <- param[4]
  rate_log <- param[5]
  
  # transform to regular parameters
  pi0 <- exp(beta0) / (1 + exp(beta0))
  pi1 <- exp(beta) / (1 + exp(beta))
  mu <- exp(mu_log)
  sigma <- exp(sigma_log)
  rate <- exp(rate_log)
  
  #combination <- delay(drug_history) + decay(drug_history)
  #combination / max(combination)
  
  # to make sure that the highest value is indeed 1
  normalizing_factor <- dnorm(mu, mu, sigma)
  
  # determine the risks given the time_since last exposure and the rate 
  # of the withdrawal model
  freq_table <- freq_table %>% mutate(
    #risk_value = exp(-rate * (unique_value - 1))
    risk_value_delayed = dnorm(unique_value, mu, sigma) / normalizing_factor, 
    risk_value_decaying = exp(-rate * (unique_value - 1)), 
    risk_value = (risk_value_delayed + risk_value_decaying) / max(risk_value_delayed + risk_value_decaying)
  )
  
  # the risk when the time_since is 0 (never exposed or currently exposed)
  # is zero
  freq_table$risk_value[1] <- 0
  
  # determine the log-likelihood given the risks
  freq_table <- freq_table %>% rowwise() %>% 
    mutate(loglikelihood = -1*(n_adr*log((pi1 - pi0)*risk_value + pi0) + 
                                 n_no_adr * log(1 - (pi1 - pi0)*risk_value - pi0)))
  
  sum(freq_table$loglikelihood)
}


#risk_model_long_term <- function(rate, delay)

#' @export
loglikelihood_long_term <- function(param,
                                    freq_table) {
  
  # extract parameters
  beta0 <- param[1]
  beta <- param[2]
  rate_log <- param[3]
  delay_log <- param[4]
  
  # transform to regular parameters
  pi0 <- exp(beta0) / (1 + exp(beta0))
  pi1 <- exp(beta) / (1 + exp(beta))
  rate <- exp(rate_log)
  delay <- exp(delay_log)
  
  # determine the risks given the time_since last exposure and the rate 
  # of the withdrawal model
  freq_table <- freq_table %>% mutate(
    #risk_value = exp(-rate * (unique_value - 1))
    risk_value = 1 / (1 + exp(-rate * (unique_value - delay)))
  )
  
  # the risk when the time_since is 0 (never exposed or currently exposed)
  # is zero
  freq_table$risk_value[1] <- 0
  
  # determine the log-likelihood given the risks
  freq_table <- freq_table %>% rowwise() %>% 
    mutate(loglikelihood = -1*(n_adr*log((pi1 - pi0)*risk_value + pi0) + 
                                 n_no_adr * log(1 - (pi1 - pi0)*risk_value - pi0)))
  
  sum(freq_table$loglikelihood)
}
