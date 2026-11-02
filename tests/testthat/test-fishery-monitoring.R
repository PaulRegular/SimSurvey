test_that("monitoring settings leave all fishing truth unchanged", {
  withr::local_seed(100)
  sim <- fishery_population(years = 2001:2003)
  set.seed(101)
  full <- fishery_run(sim, reporting_factor = 2, reporting_cv = 0.3,
                       lengths_cap = 12, ages_cap = 2, light = FALSE)$fishery
  set.seed(101)
  none <- fishery_run(sim, reporting_factor = 0.5, sampling_coverage = 0,
                       lengths_cap = 0, light = FALSE)$fishery
  for (name in c("joint_catch", "mortality_ledger", "reference", "sp_removals")) {
    expect_equal(full[[name]], none[[name]])
  }
  truth_columns <- c("fishery", "sim", "year", "operation", "cell", "n", "effort")
  expect_equal(full$setdet[, ..truth_columns], none$setdet[, ..truth_columns])
  expect_equal(nrow(none$samp), 0)
  expect_equal(none$setdet$reported_n, none$setdet$n * 0.5)
  expect_equal(none$annual$reported_cpue, none$annual$cpue * 0.5)
  expect_equal(full$setdet$n_measured, rep(12, 3))
  expect_true(all(full$samp$age_probability > 0 & full$samp$age_probability <= 1))
  expect_true(all(full$samp$operation_probability == 1))
})

test_that("reporting uses common number and weight multipliers, preserving zeros", {
  withr::local_seed(102)
  result <- fishery_run(weight = function(length) 2, reporting_factor = 3,
                         reporting_cv = 0.4)$fishery
  expect_equal(result$setdet$reported_kg, 2 * result$setdet$reported_n)
  expect_equal(result$annual$reported_biomass_cpue, 2 * result$annual$reported_cpue)
  zero <- sim_fishery(fishery_population(), 0.5,
                       function(available, operation) available$N * 0,
                       max_operations = 2, reporting_cv = 1)$fishery
  expect_equal(zero$setdet$reported_n, c(0, 0))
  expect_equal(zero$setdet$n_measured, c(0, 0))
  no_catch <- sim_fishery(fishery_population(), 0, capture_all,
                           weight = function(length) 2)$fishery
  expect_equal(no_catch$annual$kg, 0)
})

test_that("age samples are nested and capped for both strategies", {
  withr::local_seed(103)
  random <- fishery_run(lengths_cap = 10, ages_cap = 3, age_sampling = "random",
                         light = FALSE)$fishery
  expect_equal(sum(random$samp$aged), 3)
  expect_equal(random$samp$age_probability, rep(0.3, 10))
  expect_equal(random$samp$age[random$samp$aged], random$sample_truth$true_age[random$samp$aged])
  stratified <- fishery_run(lengths_cap = 30, ages_cap = 2,
                             age_sampling = "stratified")$fishery
  aged <- stratified$samp[stratified$samp$aged, ]
  expect_true(all(table(aged$length) <= 2))
  expect_true(all(stratified$samp$measured))
  expect_equal(nrow(fishery_run(lengths_cap = 0)$fishery$samp), 0)
  expect_true(all(is.na(fishery_run(ages_cap = 0)$fishery$samp$age)))
})

test_that("preferential inclusion probabilities have the specified expected coverage", {
  p <- .sampling_probabilities(c(0, 1, 10, 100), 0.5, 2)
  expect_equal(sum(p), 2)
  expect_true(all(diff(p) >= 0 & p[-1] <= 1))
  expect_equal(.sampling_probabilities(rep(0, 4), 0.3, 10), rep(0.3, 4))
  expect_equal(.sampling_probabilities(1:4, 0.4, 0), rep(0.4, 4))
  expect_equal(.sampling_probabilities(1:4, 0, 10), rep(0, 4))
  expect_equal(.sampling_probabilities(1:4, 1, 10), rep(1, 4))
})

test_that("year controls align by names and callbacks receive resolved gear", {
  withr::local_seed(104)
  sim <- fishery_population(years = 2001:2003)
  result <- sim_fishery(sim, c(`2003` = 0.3, `2001` = 0.1, `2002` = 0.2), capture_all,
                         selectivity = function(x) 1, subset_cells = cell == 1,
                         gear_amount = function(operation) operation$sim + 1,
                         duration = function(operation) operation$gear_amount * 2,
                         n_sims = 2, sampling_coverage = 0)$fishery
  expect_equal(result$annual$target, rep(c(10, 20, 30), 2))
  expect_equal(result$setdet$effort, c(rep(8, 3), rep(18, 3)))
  expect_error(sim_fishery(sim, c(a = 0.1, b = 0.2, c = 0.3), capture_all), "names")
})

test_that("length selection controls both reference targets and captured lengths", {
  withr::local_seed(105)
  result <- sim_fishery(fishery_population(), 0.5, capture_all,
                         selectivity = function(length) as.numeric(length == 20),
                         selectivity_by = "length", subset_cells = cell == 1,
                         light = FALSE)$fishery
  expect_equal(result$annual$target, 25)
  expect_equal(result$joint_catch$length, 20)
  expect_equal(result$joint_catch$age, 2)
  expect_equal(result$reference$catch_reference, c(0, 25))
})

test_that("lognormal reporting has the specified mean and CV parameterization", {
  withr::local_seed(106)
  sim <- fishery_population()
  set.seed(107)
  baseline <- fishery_run(sim, sampling_coverage = 0)$fishery
  # With constant reporting and zero coverage the final RNG state is fishing's.
  expected <- stats::rlnorm(1, log(1.5) - log(1 + 0.4^2) / 2,
                             sqrt(log(1 + 0.4^2)))
  set.seed(107)
  noisy <- fishery_run(sim, reporting_factor = 1.5, reporting_cv = 0.4,
                        sampling_coverage = 0)$fishery
  expect_equal(noisy$setdet$reporting_multiplier, expected)
  expect_equal(noisy$setdet$reported_n, baseline$setdet$n * expected)
})
