% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit-model.R
\name{fit_model}
\alias{fit_model}
\title{Fit a Model to Cohort Data}
\usage{
fit_model(
  pair,
  model = c("no-association", "current-use", "past-use", "withdrawal", "delayed",
    "decaying", "delayed+decaying", "long-term"),
  zero_patients = 0,
  zero_timepoints = 0,
  method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN", "Brent"),
  maxiter = 1000,
  parameters = list(),
  mc.cores = 1
)
}
\arguments{
\item{model}{Label for a risk model. Can be either ...}

\item{method}{Methods used by the base \code{\link{optim}}-solver}

\item{maxiter}{Maximum number iterations for the \code{\link{optim}}-solver}

\item{cohort}{A cohort dataset. See details below.}

\item{start}{Starting point for the base \code{\link{optim}}-solver
(Default: \code{c(-1,1)})}

\item{control}{List used by the base \code{\link{optim}}-solver}
}
\value{
A model fit. A list with the items
\item{\code{n_patients}}{Total number of patients in the dataset}
\item{\code{simulation_time}}{Total number of time points}
\item{\code{loglikelihood}}{The maximal log-likelihood value}
\item{\code{beta0}}{Intercept for the logistic model}
\item{\code{beta}}{Parameter for the logistic model}
\item{\code{prob_no_adr_with_drug}}{Estimated probability
of no ADR occurring}
\item{\code{prob_ADR_with_drug}}{Estimated probability
of ADR occurring}
\item{\code{convergence}}{Convergence value of \code{\link{optim}}.
\code{0} means that the algorithm converged}
}
\description{
\code{fit_model} fits a specified risk model \code{risk_model}
to cohort data.
}
\section{Patient models}{

\code{fit_model} can currently only deal with the
\code{\link{patient_model_uninformative}}. More complex
patient models, such as \code{\link{patient_model_sex}}
require the estimation of more parameters which is currently
not supported.
}

\section{Cohort data object}{

Minimal requirement is that the cohort is a list with
two matrix with the items
\itemize{
\item{\code{drug_history}  A binary matrix of size
\code{n_patients x simulation_time}  describing
the drug prescriptions.}
\item{\code{adr_history}  A binary matrix of size
\code{n_patients x simulation_time} describing
ADR histories}
}
See \code{\link{check_cohort}} and \code{\link{generate_cohort}}
for more details on the \code{cohort} object.
}

\examples{
cohort <- expard::generate_cohort(n_patients = 1000, 
                                  simulation_time = 100, 
                                  risk_model = expard::risk_model_immediate(), 
                                  verbose = TRUE, 
                                  min_chance_drug = probability_model_constant(.3), 
                                  max_chance_drug = probability_model_constant(.7), 
                                  min_chance_adr = probability_model_constant(.3),
                                  max_chance_adr = probability_model_constant(.6))
# fit the no effect model
fit_model(cohort, risk_model = expard::risk_model_no_effect()) 

# fit the true immediate effect model
fit_model(cohort, risk_model = expard::risk_model_immediate())
# note that the estimators are close to the truth (.3 and .6) 
}
\seealso{
\code{\link{check_cohort}},\code{\link{generate_cohort}}
}
