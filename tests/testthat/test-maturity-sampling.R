test_that("survey maturity protocols mask the same latent measured sample", {
  withr::local_seed(410)
  sim <- small_population()
  sim$sim_maturity <- function(length) rep(0.5, length(length))
  run <- function(protocol) sim_survey(sim, q = function(age) 1, trawl_dim = c(1, 5),
    min_sets = 2, set_den = 0, lengths_cap = 15, ages_cap = 3,
    age_sampling = "random", light = FALSE, maturity_sampling = protocol)
  set.seed(411)
  aged <- run("aged")
  rng <- .Random.seed
  set.seed(411)
  measured <- run("measured")
  expect_identical(.Random.seed, rng)
  set.seed(411)
  none <- run("none")
  expect_identical(.Random.seed, rng)
  expect_equal(aged$samp$maturity_observed, aged$samp$aged)
  expect_equal(measured$samp$maturity_observed, measured$samp$measured)
  expect_true(all(is.na(none$samp$mature)))
  expect_true(all(!none$samp$maturity_observed))
  expect_true(all(is.na(aged$samp$mature[!aged$samp$aged])))
  expect_true(all(aged$samp$mature[aged$samp$aged] %in% 0:1))
  expect_equal(aged$samp$mature[aged$samp$aged], measured$samp$mature[aged$samp$aged])
  expect_equal(aged$sample_truth$true_mature, measured$sample_truth$true_mature)
  expect_equal(aged$sample_truth$true_mature, none$sample_truth$true_mature)
  columns <- c("set", "id", "age", "length", "measured", "aged")
  expect_equal(aged$samp[, ..columns], measured$samp[, ..columns])
  expect_equal(aged$samp[, ..columns], none$samp[, ..columns])
  expect_equal(aged$full_setdet, measured$full_setdet)
  expect_equal(aged$setdet$n_maturity, aged$setdet$n_aged)
  expect_equal(measured$setdet$n_maturity, measured$setdet$n_measured)
  expect_equal(sum(aged$samp_totals$n_maturity), sum(aged$samp$aged))
  expect_equal(sum(none$samp_totals$n_maturity), 0)
  expect_null(measured$samp$true_mature)
  expect_true(all(aged$sample_truth$measured))
})

test_that("maturity is Bernoulli, not a threshold at the midpoint", {
  withr::local_seed(412)
  sim <- small_abundance()
  sample <- data.table::data.table(age = rep(5, 20), length = rep(47.81, 20),
                                    measured = TRUE, aged = TRUE)
  set.seed(413)
  expected <- stats::rbinom(20, 1, 0.5)
  set.seed(413)
  observed <- .maturity_observations(sample, sim, "measured")
  expect_equal(observed$mature, expected)
  expect_null(sample$mature)
  sim$sim_maturity <- function(length) 0
  expect_equal(.maturity_observations(sample, sim, "measured")$mature, rep(0L, 20))
  sim$sim_maturity <- function(length) 1
  expect_equal(.maturity_observations(sample, sim, "measured")$mature, rep(1L, 20))
})

test_that("age-based fishery maturity can be observed without ageing the fish", {
  withr::local_seed(414)
  sim <- fishery_population()
  sim$sim_maturity <- function(age) as.numeric(age == 2)
  sim$maturity_by <- "age"
  result <- fishery_run(sim, lengths_cap = 20, ages_cap = 0,
                         maturity_sampling = "measured", light = FALSE)$fishery
  expect_true(all(is.na(result$samp$age)))
  expect_true(all(result$samp$maturity_observed))
  expect_equal(result$samp$mature, as.integer(result$sample_truth$true_age == 2))
  expect_equal(result$setdet$n_maturity, 20)
  expect_equal(result$annual$n_maturity, 20)
  expect_equal(result$annual$n_aged, 0)
})

test_that("fishery protocols preserve reporting, all existing samples and later RNG", {
  withr::local_seed(415)
  sim <- fishery_population(years = 2001:2003)
  sim$sim_maturity <- function(length) rep(0.5, length(length))
  a <- sim_fishery(target = 0.2, catch_model = sim_capture(1),
                    reporting_cv = 0.2, max_operations = 8, lengths_cap = 10,
                    ages_cap = 2, light = FALSE)
  b <- a
  set.seed(416)
  aged <- sim_fisheries(sim, a = a, b = b)$fishery
  rng <- .Random.seed
  a$maturity_sampling <- "measured"
  set.seed(416)
  measured <- sim_fisheries(sim, a = a, b = b)$fishery
  expect_identical(.Random.seed, rng)
  a$maturity_sampling <- "none"
  set.seed(416)
  none <- sim_fisheries(sim, a = a, b = b)$fishery
  expect_identical(.Random.seed, rng)
  expect_equal(aged$joint_catch, measured$joint_catch)
  expect_equal(aged$mortality_ledger, none$mortality_ledger)
  columns <- setdiff(names(aged$setdet), "n_maturity")
  expect_equal(aged$setdet[, ..columns], measured$setdet[, ..columns])
  expect_equal(aged$setdet[, ..columns], none$setdet[, ..columns])
  columns <- setdiff(names(aged$samp), c("mature", "maturity_observed"))
  expect_equal(aged$samp[, ..columns], measured$samp[, ..columns])
  expect_equal(aged$samp[, ..columns], none$samp[, ..columns])
  expect_equal(aged$sample_truth$true_mature, measured$sample_truth$true_mature)
  expect_equal(aged$sample_truth$true_mature, none$sample_truth$true_mature)
  expect_equal(measured$samp$mature[aged$samp$maturity_observed],
               aged$samp$mature[aged$samp$maturity_observed])
  rows <- none$samp$fishery == "a"
  expect_true(all(is.na(none$samp$mature[rows])))
})

test_that("survey programs finish all existing samples before maturity", {
  withr::local_seed(417)
  sim <- small_population()
  a <- sim_survey(trawl_dim = c(1, 5), min_sets = 2, set_den = 0,
                   lengths_cap = 20, ages_cap = 2, light = FALSE)
  set.seed(418)
  aged <- sim_surveys(sim, a = a, b = a)
  rng <- .Random.seed
  a$maturity_sampling <- "measured"
  set.seed(418)
  measured <- sim_surveys(sim, a = a, b = a)
  expect_identical(.Random.seed, rng)
  for (label in c("a", "b")) {
    one <- aged$survey_programs[[label]]
    two <- measured$survey_programs[[label]]
    columns <- c("set", "id", "age", "length", "measured", "aged")
    expect_equal(one$samp[, ..columns], two$samp[, ..columns])
    expect_equal(one$sample_truth$true_mature, two$sample_truth$true_mature)
    expect_equal(one$full_setdet, two$full_setdet)
  }
  selected <- select_survey(measured, "b")
  expect_equal(selected$sample_truth, measured$survey_programs$b$sample_truth)
  expect_equal(selected$B, measured$B)
  expect_equal(selected$S, measured$S)
})

test_that("empty samples, invalid protocols and old objects are explicit", {
  withr::local_seed(419)
  empty <- small_survey(q_by = "age", lengths_cap = 0, ages_cap = 0)
  expect_equal(sum(empty$setdet$n_maturity), 0)
  expect_true(all(is.na(empty$samp$mature)))
  expect_null(empty$sample_truth)
  fish <- sim_fishery(fishery_population(), 0, capture_all)$fishery
  expect_equal(nrow(fish$samp), 0)
  expect_equal(fish$annual$n_maturity, 0)
  expect_error(small_survey(maturity_sampling = "all"), "maturity_sampling")
  expect_error(fishery_run(maturity_sampling = NA_character_), "maturity_sampling")
  legacy <- small_population()
  legacy[c("sim_weight", "sim_maturity", "maturity_by", "B", "B_at_length", "S", "S_at_length")] <- NULL
  legacy$sp_N$B <- legacy$sp_N$S <- NULL
  expect_error(sim_survey(legacy), "Maturity model missing")
  numbers <- sim_survey(legacy, maturity_sampling = "none", trawl_dim = c(1, 5), min_sets = 2, set_den = 0)
  expect_true(all(!numbers$samp$maturity_observed))
  expect_null(numbers$B)
  legacy$sim_weight <- sim_weight()
  legacy$sim_maturity <- sim_maturity()
  legacy$maturity_by <- "length"
  upgraded <- round_sim(legacy)
  expect_equal(upgraded$N, round_sim(numbers)$N)
  expect_equal(colSums(upgraded$B), colSums(upgraded$B_at_length))
})

test_that("fisheries inherit one population weight without changing number targets", {
  withr::local_seed(420)
  sim <- fishery_population()
  config <- sim_fishery(target = 0.2, catch_model = capture_all, light = FALSE)
  set.seed(421)
  first <- sim_fisheries(sim, a = config, b = config)
  sim$sim_weight <- function(length) 2
  set.seed(421)
  constant <- sim_fisheries(sim, a = config, b = config)
  expect_equal(constant$fishery$setdet$kg, constant$fishery$setdet$n * 2)
  expect_equal(first$fishery$annual$target, constant$fishery$annual$target)
  expect_equal(first$fishery$setdet$n, constant$fishery$setdet$n)
  expect_equal(constant$B, constant$N * 2)
  expect_false("weight" %in% names(formals(sim_fishery)))
  expect_false("weight" %in% names(formals(sim_fisheries)))
})
