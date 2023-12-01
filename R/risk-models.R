##############################################
# risk-models.R 
#
# Contains all the risk models. 
# They model how the risk of suffering the ADR
# changes with the history of drug history
# of the patient. 
# The risk models return values between 0 and 
# 1. Zero is minimal risk, 1 is full risk. 
##############################################

# 
risk_models <- c('no-association', 
                 'current-use', 
                 'past-use', 
                 'withdrawal', 
                 'delayed',
                 'decaying', 
                 'delayed+decaying', 
                 'long-term')


#' Risk Model 'No Association'
#' 
#' A risk model reflects how the probability of 
#' suffering the ADR changes with the drug exposures
#' of a patient. It returns \code{0} when the drug 
#' prescription history has no effect, and \code{1} 
#' when the patient is at maximal risk.  
#' In this case, the drug has no effect on the ADR
#' what so ever.
#' 
#' @return A risk model
#' @family Risk models
#' @examples 
#' drug_history <- c(1, 0, 1, 0, 0)
#' risk_model <- risk_model_no_effect() 
#' risk_model() 
#' # -> 0
#' @export
risk_model_no_association <- function() { 
  function(drug_history, ...) { 
    rep(0, length(drug_history))
  }
}

#' Risk Model 'Immediate'
#' 
#' A risk model reflects how the probability of 
#' suffering the ADR changes with the drug exposures
#' of a patient. It returns \code{0} when the drug 
#' prescription history has no effect, and \code{1} 
#' when the patient is at maximal risk.  
#' In this case, the probability of the ADR 
#' is maximal only when the patient is currently
#' exposed.
#' 
#' @return A risk model
#' @family Risk models
#' @examples 
#' drug_history <- c(1, 0, 1, 0, 0)
#' risk_model <- risk_model_immediate() 
#' risk_model(drug_history) 
#' @export
risk_model_current_use <- function() {
  function(drug_history, ...) { 
    return(drug_history)
  }
}


#' Risk Model 'Withdrawal'
#' 
#' A risk model reflects how the probability of 
#' suffering the ADR changes with the drug exposures
#' of a patient. It returns \code{0} when the drug 
#' prescription history has no effect, and \code{1} 
#' when the patient is at maximal risk.
#' In this case, once the patient is no longer exposed
#' to the drug, the probability of the ADR peaks 
#' and then dissipates exponentially. 
#' Let \eqn{\tau} be the number of time points since
#' the patient was prescribed the drug last. The return
#' value is the given by 
#' \deqn{\exp(-\gamma \cdot (\tau - 1))} 
#' where \eqn{\gamma} is \code{rate}.
#' 
#' @param rate The rate with which the risk dissipates
#' 
#' @return A risk model
#' @family Risk models
#' @examples 
#' drug_history <- c(1, 0, 1, 0, 0)
#' 
#' risk_model <- risk_model_withdrawal(rate = 1.2) 
#' risk_model(drug_history) 
#' @export
risk_model_past <- function(past) {
  
  # check correctness input
  if (past <= 0) { 
    stop("past should be > 0") 
  }
  
  function(drug_history, ...) { 
    
    simulation_time <- length(drug_history) 
    
    sapply(1:simulation_time, function(t) { 
      as.numeric( any(drug_history[max(1,t-past):t] != 0))
    })
  }
}






#' Risk Model 'Withdrawal'
#' 
#' A risk model reflects how the probability of 
#' suffering the ADR changes with the drug exposures
#' of a patient. It returns \code{0} when the drug 
#' prescription history has no effect, and \code{1} 
#' when the patient is at maximal risk.
#' In this case, once the patient is no longer exposed
#' to the drug, the probability of the ADR peaks 
#' and then dissipates exponentially. 
#' Let \eqn{\tau} be the number of time points since
#' the patient was prescribed the drug last. The return
#' value is the given by 
#' \deqn{\exp(-\gamma \cdot (\tau - 1))} 
#' where \eqn{\gamma} is \code{rate}.
#' 
#' @param rate The rate with which the risk dissipates
#' 
#' @return A risk model
#' @family Risk models
#' @examples 
#' drug_history <- c(1, 0, 1, 0, 0)
#' 
#' risk_model <- risk_model_withdrawal(rate = 1.2) 
#' risk_model(drug_history) 
#' @export
risk_model_duration <- function(duration) {
  
  # check correctness input
  if (duration <= 0) { 
    stop("past should be > 0") 
  }
  
  function(drug_history, ...) { 
    
    simulation_time <- length(drug_history) 
    
    sapply(1:simulation_time, function(t) { 
      as.numeric( sum(drug_history[1:t]) >= duration )
    })
  }
}











#' Risk Model 'Withdrawal'
#' 
#' A risk model reflects how the probability of 
#' suffering the ADR changes with the drug exposures
#' of a patient. It returns \code{0} when the drug 
#' prescription history has no effect, and \code{1} 
#' when the patient is at maximal risk.
#' In this case, once the patient is no longer exposed
#' to the drug, the probability of the ADR peaks 
#' and then dissipates exponentially. 
#' Let \eqn{\tau} be the number of time points since
#' the patient was prescribed the drug last. The return
#' value is the given by 
#' \deqn{\exp(-\gamma \cdot (\tau - 1))} 
#' where \eqn{\gamma} is \code{rate}.
#' 
#' @param rate The rate with which the risk dissipates
#' 
#' @return A risk model
#' @family Risk models
#' @examples 
#' drug_history <- c(1, 0, 1, 0, 0)
#' 
#' risk_model <- risk_model_withdrawal(rate = 1.2) 
#' risk_model(drug_history) 
#' @export
risk_model_withdrawal <- function(rate) {
  
  # check correctness input
  if (rate <= 0) { 
    stop("rate should be > 0") 
  }
  
  function(drug_history, ...) { 
    
    simulation_time <- length(drug_history)
    
    # never took the drug
    if (sum(drug_history) == 0) { 
      return(rep(0, simulation_time))
    }
    
    
    sapply(1:simulation_time, function(t) { 
      # currently exposed or did not take the drug yet
      if (drug_history[t] == 1 || sum(drug_history[1:t]) == 0) { 
        return(0)  
      } else {
        time_steps_ago = t - max(which(drug_history[1:t] == 1))
        return(exp(-rate * (time_steps_ago - 1))) 
      }
    })
  }
}





#' Risk Model 'Delayed'
#' 
#' @param rate The rate with which the risk dissipates
#' 
#' @return A risk model
#' @family Risk models
#' @examples 
#' drug_history <- c(1, 0, 1, 0, 0)
#' 
#' risk_model <- risk_model_withdrawal(rate = 1.2) 
#' risk_model(drug_history) 
#' @export
risk_model_delayed <- function(mu, sigma) {
  
  # check correctness input
  if (mu <= 0) { 
    mu <- 1e-10
    #stop("mu should be > 0") 
  }
  
  if (sigma <= 0) { 
    stop("sigma should be > 0") 
  }
  
  # to make sure that the highest value is indeed 1
  normalizing_factor <- dnorm(mu, mu, sigma)
  
  function(drug_history, ...) { 
    
    simulation_time <- length(drug_history)
    
    # never took the drug
    if (sum(drug_history) == 0) { 
      return(rep(0, simulation_time))
    }
    
    
    sapply(1:simulation_time, function(t) { 
      # currently exposed or did not take the drug yet
      if (sum(drug_history[1:t]) == 0) { 
        return(0)  
      } else {
        time_steps_ago = t - min(which(drug_history[1:t] == 1))
        return(dnorm(time_steps_ago, mu, sigma) / normalizing_factor) 
      }
    })
  }
}









#' Risk Model 'Decaying'
#' 
#' @param rate The rate with which the risk dissipates
#' 
#' @return A risk model
#' @family Risk models
#' @examples 
#' drug_history <- c(1, 0, 1, 0, 0)
#' 
#' risk_model <- risk_model_withdrawal(rate = 1.2) 
#' risk_model(drug_history) 
#' @export
risk_model_decaying <- function(rate) {
  
  # check correctness input
  if (rate <= 0) { 
    stop("rate should be > 0") 
  }
  
  #normalizing_factor <- exp(-rate * (time_steps_ago - 1))
  
  function(drug_history, ...) { 
    
    simulation_time <- length(drug_history)
    
    # never took the drug
    if (sum(drug_history) == 0) { 
      return(rep(0, simulation_time))
    }
    
    
    sapply(1:simulation_time, function(t) { 
      # did not take the drug yet
      if (sum(drug_history[1:t]) == 0) { 
        return(0)  
      } else {
        time_steps_ago = t - min(which(drug_history[1:t] == 1))
        return(exp(-rate * time_steps_ago)) 
      }
    })
  }
}



#' @export
risk_model_delayed_decaying <- function(mu, sigma, rate) {
  
  # check correctness input
  if (rate <= 0) { 
    stop("rate should be > 0") 
  }
  
  if (mu <= 0) { 
    mu <- 1e-10
    #stop("mu should be > 0")  
  } 
  
  if (sigma <= 0) { 
    stop("sigma should be > 0")
  }
  
  delay <- risk_model_delayed(mu, sigma)
  decay <- risk_model_decaying(rate)
  
  #normalizing_factor <- exp(-rate * (time_steps_ago - 1))
  
  function(drug_history, ...) { 
    if (sum(drug_history) == 0) {
      return(drug_history)
    }
    combination <- delay(drug_history) + decay(drug_history)
    combination / max(combination)
  }
}





#' Risk Model 'Long Time After'
#' 
#' A risk model reflects how the probability of 
#' suffering the ADR changes with the drug exposures
#' of a patient. It returns \code{0} when the drug 
#' prescription history has no effect, and \code{1} 
#' when the patient is at maximal risk.
#' In this case, the risk of the ADR increases
#' only after a certain moment in time before it 
#' reaches the maximal risk.
#' Let \eqn{t_0} be the time point that the drug
#' was prescribed for the first time, and \eqn{\delta}
#' be the number of time points since first description
#' that the function returns 0.5 ("half way"). The 
#' return value is 
#' \deqn{1 / (1 + \exp(-\gamma \cdot (t_0 - \delta)))} 
#' where \eqn{\gamma} is \code{rate}. Note that this 
#' is a sigmoid function.  
#' 
#' @param rate The rate with which the risk increases
#' @param delay How long it takes before the risk is .5
#' 
#' @return A risk model
#' @family Risk models
#' @examples 
#' drug_history <- c(1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1)
#' 
#' risk_model <- risk_model_long_time_after(rate = .5, delay = 9) 
#' risk_model(drug_history) 
#' @export
risk_model_long_term <- function(rate, delay) {
  function(drug_history, ...) { 
    
    simulation_time <- length(drug_history)
    
    # moment of first prescription
    sapply(1:simulation_time, function(t) {
      # if case the drug was never prescribed 
      if (sum(drug_history[1:t]) == 0) {
        return(0)
      } else {
        # moment of first prescription
        time_since_first_prescription = t - min(which(drug_history[1:t] == 1))
      
        # use a sigmoid function to determine the effect
        return(1 / (1 + exp(-rate * (time_since_first_prescription - delay))))
      }
    })
  }
}
