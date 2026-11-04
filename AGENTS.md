# Working on SimSurvey

## Purpose and architecture

SimSurvey compares sampling designs against known, spatially correlated,
age-structured populations. Keep the population process, sampling process, and
estimation of population characteristics distinguishable. Preserve the ability
to compare estimates with the corresponding truth.

The main workflow is `sim_abundance()` |> `sim_distribution()` |>
`sim_survey()` |> `run_strat()` |> `strat_error()`. Functions accept and extend a
named simulation list. Preserve existing components and their meanings when
adding a stage; do not silently replace survey outputs with another data source.

Configurable closures are intentional public extension points. Recruitment,
mortality, growth, spatial effects, catchability, and analysis functions can be
supplied by users. Prefer small, composable functions to hard-coded biological
relationships. Document the inputs and outputs required from a custom function.
See `vignettes/custom_closures.Rmd` and Regular et al. (2020),
<https://doi.org/10.1371/journal.pone.0232822>.

## R conventions

- Follow the surrounding style: two-space indentation, `<-`, double-quoted
  strings, snake_case function/argument names, and aligned wrapped arguments.
  Retain established scientific names such as `N`, `R`, `Z`, and `N0`.
- Use the native `|>` pipe; the supported R version starts at 4.1.0. Historical
  analyses may still use older APIs or `%>%`; do not copy those conventions into
  new package code or reformat unrelated files.
- Use the existing `data.table`, `sf`, and `stars` representations. Avoid adding
  another manipulation or spatial framework without a specific reason.
- Internal helpers commonly start with `.`. Use namespace-qualified calls for
  optional packages, guard them with `requireNamespace()`, and give an actionable
  error when unavailable. Keep optional features out of the required core path.
- Follow nearby `data.table` conventions, including local `NULL` declarations
  for non-standard evaluation where needed by package checks. Use
  `data.table::copy()` when new code would otherwise modify shared input tables
  by reference; check that repeated scenarios do not contaminate one another.
- Do not set hidden seeds, register persistent global options, or change the
  user's working directory in simulation functions. Clean up resources created
  by new code, including temporary files and parallel clusters, on error as well
  as success.

## Scientific and data contracts

- Preserve age/year matrix dimensions and dimnames, cell/set identifiers, and
  alignment between arrays, tables, and grids. Do not rely on incidental row
  order after joins. Keep dimensions when subsetting matrices and arrays.
- Check conservation where appropriate: distributing abundance across cells or
  converting age to length must preserve totals. Rounding changes the population
  represented by `round_sim()`; its returned summaries must match its cells.
- State units and distinguish abundance, density, catch, biomass, probabilities,
  and instantaneous rates. Grid distances/areas and fishing dimensions must use
  compatible units. Do not infer kilograms from counts without a weight model.
- Preserve the package's DFO length-group labels and consistent length groups
  between simulated truth, samples, age-length keys, and stratified analyses.
- Validate documented domains and dimensions with informative errors. Test
  mathematical behavior independently rather than copying implementation code.
- Keep scientific model changes explicit. A failing test is not permission to
  redesign mortality, movement, selectivity, or an estimator. Fix clear defects
  narrowly and record broader questions for discussion.

## Documentation and repository layout

- Edit roxygen comments in `R/` for function documentation. Include assumptions,
  argument and output meanings, units, extension contracts, and small examples.
  Regenerate the affected `man/` files and `NAMESPACE` with roxygen2; do not hand
  edit generated documentation. Avoid unrelated generator/version churn.
- Record user-visible changes and bug fixes in `NEWS.md`. Keep development
  versions distinct from releases; current work toward 0.2.0 uses 0.1.8.9000.
- Verify vignette resources inside the built source archive, not only a local
  render. Embed generated local widget resources because R cleans their support
  directories; explicitly mark web-loaded scripts external to avoid build downloads.
- Keep research analyses, exploratory scripts, and large simulations under
  `analysis/` (excluded from builds). The archived informal tests are historical
  examples, not the automated test suite. `data-raw/` contains data preparation;
  do not replace bundled data without a reproducible preparation change.
- Do not commit local R libraries, check/build results, credentials, or generated
  research outputs. Keep changes focused and preserve existing user work.

## Tests and checks

The automated suite uses testthat edition 3 in `tests/testthat/`.

- Run `testthat::test_local()` during development. Use tiny local fixtures and
  `withr::local_seed()` for stochastic tests; reset each test's RNG explicitly.
- Prefer analytical expectations, conservation, dimensions, sampling bounds,
  documented errors/warnings, and same-seed reproducibility to large snapshots
  or fragile Monte Carlo thresholds. Include zero catch and sampling boundaries
  when changing sampling code.
- Tests must not access the network, launch interactive displays, write outside
  temporary directories, or require INLA for the core simulation workflow.
  Explicitly skip optional-feature tests when their dependencies are absent.
- Test behavior changes and add a regression test for each defect fix. Exercise
  custom closures as well as defaults. Do not equate equal seeds with identical
  random results across serial and parallel implementations.
- Before considering work complete, build the source package and run
  `R CMD check --no-manual` on it. Confirm that `tests/testthat.R` is included and
  executed. The GitHub Actions matrix checks Windows, macOS, and multiple R
  versions. Report the checks actually run and any remaining warnings or skips.

## Fisheries development

Fisheries attribute existing total mortality to fishing; they do not subtract
deaths again from future pre-generated populations. Keep stock and catch tables
aggregate, with exact whole-fish mortality conservation. Complete all fishing
truth before monitoring so monitoring controls cannot alter fishing outcomes.
Keep effort units explicit and distinguish true catches from reported estimates.
Weight and maturity belong to the population. Refresh `B` and `S` from the joint
length-age key when abundance changes; do not substitute weight at mean length
or multiply separately averaged weight and maturity. `S` is all-fish expected
mature biomass, without spawning timing or sex structure. Complete every existing
sample before maturity draws, and draw latent maturity for every measured fish
before applying observation masks. Preserve the established `sim_length` component.
Batch functions take `sim` first and named configurations thereafter. Declare
all fisheries and their user-supplied allocations upfront, interleave operations
on one stock/budget, and generate survey-program truth independently before
biological subsampling.
Use only the current year's spatial population when initializing fishery stocks;
test starting ledger counts against that year's cells, not just ledger balance.
Keep `survey_programs` separate from the `surveys` scenario-design table, and
restore matching spatial availability when extracting a survey for analysis.
See `vignettes/fisheries.Rmd` for the implemented model and limitations;
`analysis/fisheries-design.Rmd` remains a historical research note.
