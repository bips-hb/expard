#' @export
#' @rdname fit_model
#' @param parameters Passed to [fit_model()]
#' @param mc.cores Number of cores used for parallelization
fit_all_models <- function(
    cohort,
    models = c(
      "no-association",
      "current-use",
      "past-use",
      "withdrawal",
      "delayed",
      "decaying",
      "delayed+decaying",
      "long-term"
    ),
    zero_patients = 0,
    zero_timepoints = 0,
    method = c(
      "L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN",
      "Brent"
    ),
    maxiter = 1000,
    parameters = list(),
    mc.cores = 1
) {

  run_parallel <- FALSE

  if (mc.cores > 1) {
    run_parallel <- TRUE
    cluster <- parallel::makeCluster(mc.cores)
    parallel::clusterEvalQ(cluster, {
      require(dplyr)
      require(expard)
    })
    parallel::parallelclusterExport(
      cluster, c("pair", "models", "method", "maxiter", "parameters"),
      envir = environment()
    )

    cat(sprintf("Running in parallel with %d cores\n", mc.cores))
  }

  fit_model_local_function <- function(model) {
    cat(sprintf("Fitting model %s...\n", model))
    fit_model(cohort, model, zero_patients, zero_timepoints, method, maxiter, parameters)
  }

  if (run_parallel) {
    # if past-use is a model, we divide it into many model fits,
    # one for each value of the past parameter. This makes it easier
    # to parallelize the entire thing
    if ("past-use" %in% models) {
      models <- models[models != "past-use"]
      simulation_time <- ncol(cohort$drug_history)
      past <- 1:(simulation_time - 1)
      models <- c(models, sprintf("past-use(%d)", past))
    }

    res_temp <- parallel::parLapply(cluster, models, function(model) {
      cat(sprintf("Fitting model %s...\n", model))
      fit_model(cohort, model, zero_patients, zero_timepoints, method, maxiter, parameters)
    }) # fit_model_local_function(model))
  } else {
    res_temp <- lapply(models, function(model) {
      cat(sprintf("Fitting model %s...\n", model))
      fit_model(cohort, model, zero_patients, zero_timepoints, method, maxiter, parameters)
    }) # fit_model_local_function(model))
  }
  if (run_parallel) {
    parallel::stopCluster(cluster)
  }

  res <- res_temp[[1]]

  for (i in 2:length(res_temp)) {
    res <- dplyr::full_join(res, res_temp[[i]])
  }

  # res$BIC <- res$n_param * log(res$n_patients * res$simulation_time) + 2*res$loglikelihood

  # determine the posterior probability for each model
  # use the parameter setting of the model with the best BIC
  min_BIC <- min(res$bestBIC)

  r <- res |>
    dplyr::group_by(.data$model) |>
    dplyr::slice_min(n = 1, .data$BIC) |>
    dplyr::mutate(delta_BIC = .data$BIC - min_BIC)

  denominator <- sum(exp(-r$delta_BIC / 2))

  r <- r |>
    dplyr::mutate(
      posterior = exp(-.data$delta_BIC / 2) / denominator
    ) |>
    dplyr::select(c("model", "posterior"))

  r |>
    dplyr::full_join(res) |>
    dplyr::arrange(dplyr::desc(posterior))
}
