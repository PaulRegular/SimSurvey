# SimSurvey 0.1.8.9000 (development)

* Add fishery `target` modes for shares of Z, fully selected F (`F_full`), catch
  numbers (`catch_n`), and catch kilograms (`catch_kg`). Retain `fishery_share`
  as a deprecated alias and preserve positional calls and saved configurations.
* Validate starting-population limits, mortality-compatible targets and combined
  reference F before random stock allocation, with fishery/year/age diagnostics.
  Cap kilogram catches as whole-fish aggregate removals and record target units,
  requested settings, reference catches, achieved catches and shortfalls.
* Update fishery plots and the workflow vignette for mixed target types,
  reference F and kilogram quotas, including the mathematical conversion.
* Bound jointly rounded reference-number targets by the shared integer death
  count, preventing floating-point round trips at half-integer ties from
  creating an unsupported extra fish.

* Publish development documentation under `/dev/`, preserving the released site
  at the root. Build the site on pull requests before deployment from `main`.
  Use stable appendix anchors in both the vignette and website article.
* Shorten the fisheries vignette, replace verification tables with explanations,
  and show a function-based workflow. Extend examples to ten years, improve
  imperfect harvester knowledge, and compare random with targeted fishing using
  series divided by their overall temporal means.
* Specify the survey-map polygon trace type explicitly, removing Plotly's
  inferred-scatter message without suppressing messages from plotting.
* Base the new-workflow vignette on the documented five-age, 20-km survey example,
  with default biology, abundance, spatial processes, and survey allocation.
* Replace the vignette's fixed-core CPUE example with default shifting spatial
  distributions and imperfect current-year knowledge. Explain how changing
  concentration can steepen or flatten the CPUE-abundance relationship.
* Clarify that fishery knowledge uses current-year density and fresh annual map
  errors, and verify that fishing locations follow a shifted population.
* Use one controlled five-year population throughout the fisheries vignette,
  increase survey coverage to eight sets, match CPUE markers to their lines,
  and size widgets without deprecated Plotly layout arguments.
* Increase spacing in fishery and survey map/composition graphics and keep
  zero-catch survey locations visible, including when marker counts are constant.
* Prefer `total_deaths` in fishery reference diagnostics, clarifying that it is
  the shared age/year whole-fish total across fishing and other mortality.
  Retain `death_budget` as an identical compatibility alias.
* Embed local widget support resources in the workflow vignette so its packaged
  HTML remains portable, while keeping one Plotly web bundle external and
  avoiding plotting-asset downloads during builds.
* Revise the multi-fishery vignette around user-supplied allocations, explicit
  shortfalls/F diagnostics, survey catch kilograms, and package Plotly helpers.
  Add a controlled five-year preferential-fishing CPUE demonstration and its
  spatial-proportionality equations; separate backwards compatibility guidance.
* Add Plotly fishery maps/compositions and annual trends with `plot_fishery()`
  and `plot_fishery_trend()`. Extend `plot_trend(mat = ...)` to biomass and mature
  biomass, with matching units and preserved positional calls.
* Add deterministic individual survey catch weights and catch kilograms per set,
  age/set, and replicate/year, including batch, parallel, and scenario workflows.
  Weights use all caught fish, preserve RNG state, and remain missing when an
  older number-only simulation has no weight model.
* Retain realized catches, removal fractions, and Baranov-equivalent F at age in
  fishery reference diagnostics, alongside original N, Z, and integer budgets.
  Clarify capture q units and the allocation/shortfall interpretation.
* Correct fishery stock initialization to use only the current year's spatial
  population. Previous multi-year fishery results pooled years and must be rerun
  from their saved populations; population trajectories need not be regenerated.
* Expand "Simulating multiple fisheries and surveys" with default-based workflows,
  controlled targeting and monitoring comparisons, maturity observations,
  independent surveys, and a mathematical appendix explaining the extensions.
* Add population-level weight and maturity relationships, `sim_weight()` and
  `sim_maturity()`, and a `sim_growth()` alias for `sim_vonB()`. Derive expected
  biomass (`B`, `B_at_length`) and mature biomass (`S`, `S_at_length`) from the
  joint length-age key; refresh mass summaries in spatial and rounded populations.
* Move weight configuration from the unreleased fishery interfaces to
  `sim_abundance()`. Every fishery now shares the population weight relationship.
* Prefer `x50` to `x0` in `sim_logistic()`, retaining positional calls and a
  deprecated named `x0` alias without changing the default curve.
* Add `maturity_sampling = "aged"`, `"measured"`, or `"none"` to fisheries and
  surveys, including batch, parallel, scenario, and resume workflows. Maturity
  observations and counts mask shared latent states after all existing sampling,
  preserving catches, reporting, existing samples and overlapping maturity states.
* Add `sim_fisheries(sim, ...)` for interleaved fisheries sharing one stock and
  mortality budget, with joint target rounding and per-fishery effort/CPUE.
  Support unnamed mortality matrices returned by existing custom closures, while
  rejecting conflicting age/year labels.
* Add `sim_surveys(sim, ...)` for independent programs observing one annual
  snapshot, and `select_survey()` for existing analyses and plots. Preserve the
  scenario-design table used by `test_surveys()`. Batch configurations are named
  and canonicalized; common replicate counts align across fisheries and surveys.
* Add aggregate fishery simulation with explicit capture models, age/length gear
  selectivity, conditional mortality budgets, spatial depletion, effort and CPUE.
  Integrated reporting and nested biological sampling retain their inclusion
  probabilities and do not change fishing truth.
* Add length-based survey catchability with `q_by = "length"`, preserving the
  historical age-based default. Round and simulate on copies of caller tables.
* Preserve dimensions in single-age length keys, spatial distributions, and
  rounded single-age or single-year populations.
* Convert legacy `inla.mesh` inputs automatically for current INLA/fmesher in
  standard and barrier SPDE simulations, preserving the stored geometry and
  barrier triangle assignments. Use the current fmesher constructor in
  `make_mesh()` when available, avoiding INLA's constructor deprecation warning.
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
* Preserve the retained survey's `I_at_length` truth when running or resuming
  survey scenarios.
* Explicitly assign the documented kilometre-based projected CRS in
  `make_grid()`; the data-frame conversion could leave the grid without a CRS.

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
