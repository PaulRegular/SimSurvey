# SimSurvey 0.1.8.9000 (development)

* Add repository guidance and an automated testthat suite covering the core
  population, spatial, survey, and analysis workflow, including parallel runs
  and export/resume. Archive the historical informal tests under `analysis/`.
* Include tests in source-package checks and enable the existing CI matrix on
  the `fisheries-0.2.0` development branch.
* Reject descending, duplicated, empty, non-numeric, and non-finite age/year
  sequences in `sim_abundance()` with the documented sequence error.
* Apply the current innovation scale throughout `sim_ays_covar()` when `sd`,
  `phi_age`, or `phi_year` varies. Cache only the spatial factor, and clarify that
  `sd` is a standard deviation rather than a variance.
* Fix `strat_means()` for a single survey group, including replacement of
  undefined confidence-limit values with `NA`.
* Allow `group_lengths()` to return an empty vector, so surveys with no catch
  or no measured lengths retain their set records without failing.

# SimSurvey 0.1.1

* Initial release documented in Regular et al. (2020) <doi.org/10.1371/journal.pone.0232822>

# SimSurvey 0.1.2

* Initial release to CRAN

# SimSurvey 0.1.3

* Simplify make_grid function to improve the splits of the divisions and strata
* CRAN request: change doi in DESCRIPTION to preferred format

# SimSurvey 0.1.4

* Add a "bezier" method to make_grid
* Allow a vector of age-specific parameters to be supplied to sim_parabola plus add some options for defining a more asymmetric parabola
* Speed up sim_sets and ensure sample call passes new error traps
* Fix bug in plotting scripts; plotly returns an error if supplied xtabs class data
* Improve vis_fit function

# SimSurvey 0.1.5

* Add informative error if NaN values are generated using sim_distribution
* Fix showscale problem in sim_distribution; was supplying a logical vector when length should have been 1
* Allow for more flexible simulation of sets by adding an argument to sim_sets called subset_cells and a custom_sets argument to sim_survey
* Add number of fish available to the survey (I = N * q) to the sp_N object when running sim_survey

# SimSurvey 0.1.6

* Switch to dependence on the sf and stars package rather than the sp and raster package given the evolution of R spatial
* Switch to new INLAspacetime functions for simulating from a barrier mesh

# SimSurvey 0.1.7

* Minor improvements to documentation
* Dropped dependence on **magrittr** pipe (`%>%`); now using native R pipe (`|>`)

# SimSurvey 0.1.8

* Add `Matrix` to `Imports` so lazy-loaded mesh data (`survey_mesh`, `survey_lite_mesh`) pass checks. No user-visible changes; API and results unchanged.
