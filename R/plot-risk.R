#' Plot Risk
#'
#' Plots the risk function's behavior for
#' a given drug prescription history.
#' The \eqn{x}-axis denote the
#' drug prescriptions and the `y`-axis represent the
#' risk level, i.e., `1` is maximal risk, `0`
#' minimal risk.
#'
#' @param drug_history Binary vector denoting the
#'                    drug prescription history.
#'              (Default: `0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0`)
#' @param risk_model One of the risk models
#' @param simulation_time If not `NULL`, total number of time points, i.e.,
#'                   `0 0 0 0 1 1 1 1 1 1 0 0 0 0 ...`. (Default: `NULL`)
#' @param title Title of the plot
#' @param ylim The limits of the y-axis: `c(ymin, ymax)`
#'                (Default: `c(0,1)`)
#' @param shaded_area If `TRUE`, periods in which drugs were
#'             prescribed are shaded in a gray tone. Otherwise,
#'             the drug prescription periods are denoted by
#'             vertical dashed red lines (Default: `TRUE`)
#' @param fill Color used for the shaded area (Default: `black`)
#' @param alpha Alpha value for shaded area (Default: `.3`)
#'
#' @return A [ggplot2] plot
#' @examples
#' drug_history <- c(rep(0, 4), rep(1, 6), rep(0, 10))
#' risk_model <- risk_model_withdrawal(rate = 3)
#'
#' # create the plot
#' p <- plot_risk(drug_history,
#'   risk_model,
#'   title = "Direct effect model"
#' )
#' p
#' @export
plot_risk <- function(
    drug_history = c(rep(0, 4), rep(1, 6), rep(0, 10)),
    risk_model = risk_model_decaying(3),
    simulation_time = NULL,
    title = "",
    ylim = c(0, 1),
    shaded_area = TRUE,
    fill = "black",
    alpha = 0.3
) {

  if (!is.null(simulation_time)) {
    if (simulation_time < 10) {
      stop("simulation time must be at least 10")
    }
    drug_history <- c(rep(0, 4), rep(1, 6), rep(0, simulation_time - 10))
  }

  # determine the risks given the drug prescription history
  # and the risk model given by risk_model
  risks <- risk_model(drug_history)

  # create a dataset with the time points, the drug prescriptions
  # and the risks
  xdf <- data.frame(
    t = 1:length(drug_history),
    drug = drug_history,
    risk = risks
  )

  # determine the time points when the drug prescriptions changes from
  # "not prescribed" to "prescribed" and the other way around
  changes <- which(diff(drug_history) != 0)
  # if drug is prescribed at the first time point, add it to changes
  if (drug_history[1] == 1) {
    changes <- c(0, changes)
  }
  # if drug is prescribed at last time point, add it to changes
  if (utils::tail(drug_history, 1) == 1) {
    changes <- c(changes, length(drug_history))
  }

  # create a data frame for the change points on the x-axis. Note
  # that we set it off by + 0.5 so that it falls in the middle
  change_points <- data.frame(changes = changes + 0.5)

  # plot the risk over time
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = xdf, mapping = ggplot2::aes(x = t, y = risks)) +
    ggplot2::scale_x_continuous(limits = c(0.45, length(drug_history) + .55), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(.01, .01)) +
    ggplot2::labs(
      title = title,
      x = "time (t)",
      y = "Risk Level"
    )

  # add a shaded area when the patient is prescribed the drug
  if (shaded_area) {
    # number of periods a person was prescribed the drug
    n_drug_periods <- nrow(change_points) / 2

    if (n_drug_periods > 0) {
      temp <- data.frame(
        xmin = change_points[seq(1, nrow(change_points), by = 2), 1],
        xmax = change_points[seq(2, nrow(change_points), by = 2), 1]
      )
      p <- p +
        ggplot2::geom_rect(
          data = temp,
          ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -Inf, ymax = Inf),
          alpha = alpha,
          fill = fill
        )
    }
  } else {
    p <- p + ggplot2::geom_vline(
      data = change_points,
      ggplot2::aes(xintercept = changes),
      linetype = "dashed",
      color = "red"
    )
  }

  p <- p +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")
    )

  return(p)
}
