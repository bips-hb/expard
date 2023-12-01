# `expard`: Simulating and Fitting Exposure Models
-----------------------------------------------

`expard` constitutes an `R` package encompassing two primary functionalities:

1. It serves as a simulator for electronic healthcare data, enabling the simulation of extensive patient populations observed across multiple time points. Leveraging `R`'s functional programming features, the package facilitates the straightforward specification of intricate relationships between drug exposures and adverse drug reactions (ADRs). These relationships may be contingent on various patient attributes such as sex, age, region, among others.

2. The package includes a methodology designed to fit eight exposure models outlined in the associated paper to electronic healthcare data.

## Usage

### Simulator

The main function of the package is `generate_cohort` (see `?generate_cohort`). For example, 

```R 
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
``` 
simulates a cohort, where 

* `n_patients` is the number of patients simulated in the cohort 

* `simulation_time` is the total number of time steps

* `n_drug_ADR_pairs` is the number of drug-ADR pairs simulated

* `risk_model` Vector with risk models. Each risk model is given as a string, e.g., `"risk_model_current_use()"`. The vector must have the length `n_drug_ADR_pairs` 

* `min_chance_drug` Vector with the probabilities of the drug being prescribed
            when the drug history and the ADR history have no effect. 
            Must have a length of `n_drug_ADR_pairs`. 

* `avg_duration` Average number of time points a patient is exposed
                     once exposed to the drug (Determines `max_chance`)
                     Must have a length of `n_drug_ADR_pairs`. 

* `max_chance_drug` The probability of the ADR when 
            the drug history and the ADR history have the highest 
            possible effect (Default: `NULL`)
            Must have a length of `n_drug_ADR_pairs`. 

* `guaranteed_exposed` If `TRUE`, the patient is exposed to the drug
            at least once
            Must have a length of `n_drug_ADR_pairs`. 

* `min_chance` The probability of the ADR when 
            the drug history has no effect
            Must have a length of `n_drug_ADR_pairs`. 

* `max_chance` The probability of the ADR when 
            the drug history has the highest possible effect
            Must have a length of `n_drug_ADR_pairs`.

### Exposure Model Fit

One can fit the eight exposure models by running `fit_all_models.R`. For more 
information, see the comments to `fit_model`.

## Acknowledgements

We gratefully acknowledge the financial support from the innovation fund (“Innovationsfonds”) of the Federal Joint Committee in Germany (grant number: 01VSF16020).

## Contact

Louis Dijkstra\
Leibniz Institute for Prevention Research & Epidemiology  
E-mail: dijkstra (at) leibniz-bips.de