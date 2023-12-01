#' Plot the Fit 
#' 
#' Plots the BIC for the different risk models that were fitted by the 
#' function \code{\link{fit_all_models}}. 
#' 
#' @param fit The fit 
#' @param xlab The x-label
#' @param past_values If not \code{NULL}, a select number of values for past are 
#'             used for the plot
#' 
#' @return ggplot 
#' @export
plot_fit <- function(fit, 
                     x_label = "model", 
                     title = "", 
                     y_range = NULL, 
                     past_values = NULL) {
  
  
  # if one run the fit_all_models in parallel, process the data 
  if('past-use(1)' %in% fit$model) {
    simulation_time <- fit$simulation_time[1]
    
    past_models <- sprintf('past-use(%d)', 1:(simulation_time-1))
    
    fit$model <- sapply(1:nrow(fit), function(i) {
      if(grepl('past-use', fit$model[i], fixed = TRUE)) {
        return('past-use')
      } else {
        fit$model[i]
      }
    })
  }
  
  # Use only the 'past-use' models for which the past parameter falls in the given range
  if (!is.null(past_values)) {
    fit <- fit %>% filter(model != 'past-use' | past %in% past_values)
  }

  # get the best BIC fit for each model
  best_fit <- fit %>% group_by(model) %>% 
    filter(BIC == min(BIC))  %>% 
    arrange(past) %>% 
    filter(row_number() == 1) %>% 
    arrange(BIC)

  # get the overall minimum and maximum BIC value
  min_BIC <- min(best_fit$BIC)
  max_BIC <- max(best_fit$BIC)
  
  if (is.null(y_range)) {
    y_range <- c(min_BIC,max_BIC)
  }

  old_model_label <- c("past-use", "current-use", "no-association")
  new_model_label <- c("past use", "current use", "no association")
  
  # Use sapply to replace each occurrence in the input_vector
  best_fit$model <- sapply(best_fit$model, function(x) {
    for (i in seq_along(old_model_label)) {
      x <- gsub(old_model_label[i], new_model_label[i], x)
    }
    return(x)
  })
  
  # plot just the best fit
  ggplot(best_fit) +
    geom_bar(aes(x = reorder(model, BIC), y = BIC), stat="identity") +
    coord_cartesian(ylim=y_range) + 
    ggtitle(title) + 
    #scale_y_continuous(expand = expansion(mult = c(0.1, .1))) + #expand = c(0, 100)) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlab(x_label) 
  
  # plot just the best fit
  # ggplot(best_fit) +
  #   geom_bar(aes(x = reorder(model, BIC), y = BIC), stat="identity") +
  #   coord_cartesian(ylim=y_range) + 
  #   scale_y_continuous(expand = expansion(mult = c(0.1, .1))) + #expand = c(0, 100)) + 
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #   xlab("model") 
}
