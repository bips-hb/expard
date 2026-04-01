# expard 1.1.0

* `create2x2table()` drug-era method has been rewritten to correctly track state across time points, fixing a bug where era transitions and ADR occurrences were not properly recorded (#1).
* `create2x2table()` and `create2x2tables()` examples have been corrected to use the appropriate input types and function names (#1).
* `fit_model()` now correctly references the `cohort` parameter throughout; several internal code paths were still using a previous parameter name, causing runtime errors for most model types (#1).
* `fit_model()` examples have been corrected to use valid function and parameter names (#1).
* `fit_all_models()` fixed a typo (`parallelclusterExport` -> `clusterExport`) that prevented parallel execution from working (#1).
* Risk model examples have been corrected to use the actual exported function names (#1).
* Fixed roxygen2 documentation to remove references to non-existent functions and resolve all R CMD check warnings and notes (#1).
