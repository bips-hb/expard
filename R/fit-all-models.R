#' @export
fit_all_models <- function(pair,
                           models = c(
                             'no-association',
                             'current-use',
                             'past-use',
                             'withdrawal',
                             'delayed',
                             'decaying',
                             'delayed+decaying',
                             'long-term'
                           ),
                           zero_patients = 0, 
                           zero_timepoints = 0, 
                           method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN",
                                      "Brent"),
                           maxiter = 1000, 
                           parameters = list(),
                           mc.cores = 1) { 
  
  run_parallel <- FALSE
  if (mc.cores > 1) {
    run_parallel <- TRUE
    cluster <- makeCluster(mc.cores)
    clusterEvalQ(cluster, {library(dplyr); library(expard)})
    clusterExport(cluster, c("pair", "models", "method", "maxiter", "parameters"), envir = environment())
    
    cat(sprintf("Running in parallel with %d cores\n", mc.cores))
  }
  
  fit_model_local_function <- function(model) {
    cat(sprintf("Fitting model %s...\n", model))
    fit_model(pair, model, zero_patients, zero_timepoints, method, maxiter, parameters) 
  }
  
  if (run_parallel) {
    
    # if past-use is a model, we divide it into many model fits, 
    # one for each value of the past parameter. This makes it easier
    # to parallelize the entire thing
    if ('past-use' %in% models) {
      models <- models[models != 'past-use']
      simulation_time <- ncol(pair$drug_history) 
      past <- 1:(simulation_time - 1)
      models <- c(models, sprintf("past-use(%d)", past))
    }
    
    res_temp <- parLapply(cluster, models, function(model) { 
      cat(sprintf("Fitting model %s...\n", model))
      expard::fit_model(pair, model, zero_patients, zero_timepoints, method, maxiter, parameters) 
      })#fit_model_local_function(model))
  } else {
    res_temp <- lapply(models, function(model) { 
      cat(sprintf("Fitting model %s...\n", model))
      expard::fit_model(pair, model, zero_patients, zero_timepoints, method, maxiter, parameters) 
    })#fit_model_local_function(model))
  }
  if (run_parallel) {
    stopCluster(cluster)
  }
  
  res <- res_temp[[1]]
  
  for(i in 2:length(res_temp)) { 
    res <- full_join(res, res_temp[[i]])  
  }
  
  #res$BIC <- res$n_param * log(res$n_patients * res$simulation_time) + 2*res$loglikelihood
  
  # determine the posterior probability for each model 
  # use the parameter setting of the model with the best BIC 
  min_BIC <- min(res$bestBIC)
  
  r <-  res %>% 
    group_by(model) %>% 
    slice_min(n = 1, BIC) %>% 
    mutate(delta_BIC = BIC - min_BIC)
  
  denominator <- sum(exp(-r$delta_BIC/2))
  r <- r %>% mutate(
    posterior = exp(-delta_BIC / 2) / denominator
  ) 
  
  r <- r %>% dplyr::select(model, posterior)
  
  res <- full_join(r, res)
  
  return(res %>% dplyr::arrange(-posterior))
}