% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate-cohort.R
\name{generate_cohort}
\alias{generate_cohort}
\title{Generate a Cohort}
\usage{
generate_cohort(
  n_patients = 100,
  simulation_time = 100,
  n_drug_ADR_pairs = 50,
  risk_model = rep("risk_model_current_use()", n_drug_ADR_pairs),
  min_chance_drug = rep(0.1, n_drug_ADR_pairs),
  avg_duration = rep(5, n_drug_ADR_pairs),
  max_chance_drug = rep(NULL, n_drug_ADR_pairs),
  prob_guaranteed_exposed = rep(1, n_drug_ADR_pairs),
  min_chance = rep(0.1, n_drug_ADR_pairs),
  max_chance = rep(0.4, n_drug_ADR_pairs),
  verbose = FALSE
)
}
\arguments{
\item{n_patients}{Number of patients}

\item{simulation_time}{The total number of time steps}

\item{n_drug_ADR_pairs}{Number of drug-ADR pairs simulated}

\item{risk_model}{Vector with risk models. Each risk model is given as
a string, e.g., \code{"risk_model_current_use()"}. The
vector must have the length \code{n_drug_ADR_pairs}#'}

\item{min_chance_drug}{Vector with the probabilities of the drug being prescribed
when the drug history and the ADR history have no effect.
Must have a length of \code{n_drug_ADR_pairs}.}

\item{avg_duration}{Average number of time points a patient is exposed
once exposed to the drug (Determines \code{max_chance})
Must have a length of \code{n_drug_ADR_pairs}.}

\item{max_chance_drug}{The probability of the ADR when
the drug history and the ADR history have the highest
possible effect (Default: \code{NULL})
Must have a length of \code{n_drug_ADR_pairs}.}

\item{min_chance}{The probability of the ADR when
the drug history has no effect
Must have a length of \code{n_drug_ADR_pairs}.}

\item{max_chance}{The probability of the ADR when
the drug history has the highest possible effect
Must have a length of \code{n_drug_ADR_pairs}.}

\item{verbose}{Show progress bar (Default: \code{FALSE})}
}
\value{
A \code{cohort} object; A list of length \code{n_drug_ADR_pairs}.
Each entry contains two binary matrices; one for the drug histories
of all patients and one for the ADR histories of all patients. Each
row is a patient
}
\description{
Generates a cohort which is a collection of drug-ADR pair.
It is a wrapper function for the \code{\link{generate_patient}}
function.
}
\seealso{
\code{\link{generate_patient}}
}
