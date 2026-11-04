test_that("survey kilograms sum all caught fish before biological subsampling", {
  withr::local_seed(721)
  sim <- small_population()
  sim$sim_weight <- function(length) (length / 10)^3
  result <- sim_survey(sim, q = function(age) 1, trawl_dim = c(1, 5),
    min_sets = 2, set_den = 0, lengths_cap = 2, ages_cap = 1, light = FALSE)
  expect_equal(result$samp$weight, (result$samp$length / 10)^3)
  expected <- result$samp[, list(kg = sum(weight)), by = "set"]
  expect_equal(result$setdet$kg[match(expected$set, result$setdet$set)], expected$kg)
  expect_equal(sum(result$setdet$kg), sum(result$samp$weight))
  expect_equal(sum(result$full_setdet$kg), sum(result$samp$weight))
  expect_equal(sum(result$samp_totals$kg_caught), sum(result$samp$weight))
  expect_gt(sum(result$setdet$kg), sum(result$samp$weight[result$samp$measured]))
  expect_equal(result$sample_truth$weight, result$samp$weight[result$samp$measured])
})

test_that("survey weights do not change RNG, samples, or program isolation", {
  withr::local_seed(722)
  sim <- small_population()
  saved <- data.table::copy(sim)
  spec <- sim_survey(q = function(age) 1, trawl_dim = c(1, 5), min_sets = 2,
                      set_den = 0, lengths_cap = 20, ages_cap = 2, light = FALSE)
  set.seed(723)
  original <- sim_surveys(sim, a = spec, b = spec)
  rng <- .Random.seed
  constant <- data.table::copy(sim)
  constant$sim_weight <- function(length) 2
  set.seed(723)
  weighted <- sim_surveys(constant, a = spec, b = spec)
  expect_identical(.Random.seed, rng)
  for (label in c("a", "b")) {
    one <- original$survey_programs[[label]]
    two <- weighted$survey_programs[[label]]
    fields <- c("set", "id", "length", "age", "measured", "aged", "mature", "maturity_observed")
    expect_equal(one$samp[, ..fields], two$samp[, ..fields])
    expect_equal(one$setdet$n, two$setdet$n)
    expect_equal(two$setdet$kg, 2 * two$setdet$n)
    expect_equal(two$samp_totals$kg_caught, 2 * two$samp_totals$n_caught)
  }
  view <- select_survey(weighted, "a")
  expect_equal(view$setdet$kg, weighted$survey_programs$a$setdet$kg)
  expect_equal(view$sample_truth$weight, rep(2, nrow(view$sample_truth)))
  expect_equal(sim, saved)
})

test_that("empty catches and legacy number-only surveys have explicit weights", {
  withr::local_seed(724)
  sim <- small_population()
  legacy <- data.table::copy(sim)
  legacy[c("sim_weight", "sim_maturity", "maturity_by", "B", "S", "B_at_length", "S_at_length")] <- NULL
  legacy$sp_N$B <- legacy$sp_N$S <- NULL
  result <- sim_survey(legacy, maturity_sampling = "none", q = function(age) 1,
    trawl_dim = c(1, 5), min_sets = 2, set_den = 0, light = FALSE)
  expect_true(all(is.na(result$samp$weight)))
  expect_true(all(is.na(result$setdet$kg[result$setdet$n > 0])))
  expect_true(all(result$full_setdet$kg[result$full_setdet$n == 0] == 0))
  expect_true(all(is.na(result$samp_totals$kg_caught)))
  for (population in list(sim, legacy)) {
    empty <- sim_survey(population, maturity_sampling = "none", q = function(age) 0,
      trawl_dim = c(1, 5), min_sets = 2, set_den = 0)
    expect_equal(empty$samp$weight, numeric())
    expect_equal(empty$setdet$kg, rep(0, nrow(empty$setdet)))
    expect_equal(empty$samp_totals$kg_caught, rep(0, nrow(empty$samp_totals)))
  }
})
