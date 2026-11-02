test_that("a single fishery batch matches standalone execution", {
  withr::local_seed(301)
  sim <- fishery_population()
  config <- sim_fishery(fishery_share = 0.5, catch_model = sim_capture(2),
                         selectivity = function(x) 1, subset_cells = cell == 1,
                         max_operations = 20, light = FALSE)
  set.seed(302)
  standalone <- sim_fishery(sim, 0.5, config$catch_model, selectivity = config$selectivity,
                             subset_cells = cell == 1, max_operations = 20, light = FALSE)
  rng <- .Random.seed
  set.seed(302)
  batch <- sim_fisheries(sim, fishery = config)
  expect_equal(batch, standalone)
  expect_identical(.Random.seed, rng)
})

test_that("fisheries share mortality without double removal", {
  withr::local_seed(303)
  sim <- fishery_population()
  # Different ages may each receive 80% of their mortality: scalar shares need
  # not sum below one when their selectivities do not overlap.
  young <- sim_fishery(fishery_share = 0.8, catch_model = capture_all,
                        selectivity = function(age) as.numeric(age == 1),
                        subset_cells = cell == 1, light = FALSE)
  old <- sim_fishery(fishery_share = 0.8, catch_model = capture_all,
                      selectivity = function(age) as.numeric(age == 2),
                      subset_cells = cell == 1, light = FALSE)
  result <- sim_fisheries(sim, young = young, old = old)$fishery
  expect_equal(result$annual$target, c(40, 40))
  expect_equal(result$annual$n, c(40, 40))
  ledger <- result$mortality_ledger
  expect_equal(sum(ledger$starting), 200)
  expect_equal(sum(ledger$fishery_deaths), 80)
  expect_equal(sum(ledger$other_deaths), 20)
  expect_equal(sum(ledger$survivors), 100)
  expect_equal(ledger$starting, ledger$fishery_deaths + ledger$other_deaths + ledger$survivors)
  expect_equal(sum(result$sp_removals$n), 80)
  incompatible <- young
  incompatible$selectivity <- function(age) 1
  expect_error(sim_fisheries(sim, young = incompatible, old = old), "Combined fishery shares")
  # Validate every year, rather than only the first one.
  sim <- fishery_population(years = 2001:2003)
  incompatible$fishery_share <- c(0.1, 0.1, 0.9)
  expect_error(sim_fisheries(sim, first = incompatible, second = old), "mortality budget")
})

test_that("registration order does not change interleaved fishing or monitoring", {
  withr::local_seed(304)
  sim <- fishery_population()
  capture <- function(available, operation) .thin_counts(available$N, 2)
  first <- sim_fishery(fishery_share = 0.4, catch_model = capture,
                        selectivity = function(x) 1, subset_cells = cell == 1,
                        gear_unit = "hook", time_unit = "hour", light = FALSE)
  second <- first
  second$gear_amount <- 2
  second$gear_unit <- "net"
  set.seed(305)
  a <- sim_fisheries(sim, first = first, second = second)$fishery
  rng <- .Random.seed
  set.seed(305)
  b <- sim_fisheries(sim, second = second, first = first)$fishery
  expect_equal(a, b)
  expect_identical(.Random.seed, rng)
  sequence <- a$setdet$fishery
  expect_lt(min(which(sequence == "first")), max(which(sequence == "second")))
  expect_lt(min(which(sequence == "second")), max(which(sequence == "first")))
  for (label in c("first", "second")) {
    rows <- a$setdet[a$setdet$fishery == label, ]
    expect_equal(rows$fishery_operation, seq_len(nrow(rows)))
    annual <- a$annual[a$annual$fishery == label, ]
    expect_equal(annual$cpue, sum(rows$n) / sum(rows$effort))
  }
  expect_equal(a$setdet$gear_unit, ifelse(sequence == "first", "hook", "net"))
})

test_that("joint target rounding respects small overall mortality budgets", {
  withr::local_seed(306)
  sim <- fishery_population(3)
  config <- sim_fishery(fishery_share = 0.5, catch_model = capture_all,
                         selectivity = function(x) 1, subset_cells = cell == 1)
  result <- sim_fisheries(sim, first = config, second = config)$fishery
  # Each reference target is .75; joint rounding allocates one fish to each.
  expect_equal(result$annual$reference_target, c(0.75, 0.75))
  expect_equal(result$annual$target, c(1, 1))
  expect_equal(sum(result$annual$n), 2)
  expect_equal(sum(result$mortality_ledger$fishery_deaths), 2)
  expect_equal(sum(result$mortality_ledger$other_deaths), 0)
})

test_that("competition can exhaust selected stocks before reference targets", {
  withr::local_seed(307)
  sim <- fishery_population()
  # Fishery A is allowed at both ages but captures only age 1. A and B both
  # require that age, so their combined total target cannot be achieved.
  a <- sim_fishery(fishery_share = 0.5,
                    catch_model = function(available, operation) available$N * (available$age == 1),
                    selectivity = function(age) 1, subset_cells = cell == 1,
                    max_operations = 5)
  b <- sim_fishery(fishery_share = 0.5, catch_model = capture_all,
                    selectivity = function(age) as.numeric(age == 1), subset_cells = cell == 1)
  result <- sim_fisheries(sim, a = a, b = b)$fishery
  expect_equal(sum(result$annual$target), 75)
  expect_lte(sum(result$annual$n), 50)
  expect_gte(sum(result$annual$shortfall), 25)
  expect_true(all(result$mortality_ledger$survivors >= 0))
})

test_that("fishery batch configuration and replicate counts are validated", {
  withr::local_seed(308)
  sim <- fishery_population()
  config <- sim_fishery(fishery_share = 0.2, catch_model = capture_all)
  expect_error(sim_fisheries(sim), "nonempty")
  expect_error(sim_fisheries(sim, config), "nonempty")
  expect_error(sim_fisheries(sim, a = config, a = config), "unique")
  expect_error(sim_fisheries(sim, a = sim_survey()), "sim_fishery_spec")
  expect_error(sim_fisheries(sim, a = sim_fishery(fishery_share = 0.2,
                catch_model = capture_all, n_sims = 2)), "common n_sims")
  expect_error(sim_fisheries(sim, a = sim_fishery(fishery_share = 0.2,
                catch_model = capture_all, weight = function(length) 1)), "unused argument")
  result <- sim_fisheries(sim, a = config, n_sims = 2)
  expect_identical(result$replicates, 1:2)
  expect_error(sim_fisheries(result, a = config, n_sims = 1), "existing batch")
  inherited <- sim_fisheries(result, a = config)
  expect_equal(sort(unique(inherited$fishery$annual$sim)), 1:2)
})

test_that("monitoring changes in one fishery leave all shared truth unchanged", {
  withr::local_seed(309)
  sim <- fishery_population(years = 2001:2003)
  a <- sim_fishery(fishery_share = 0.2, catch_model = sim_capture(1), max_operations = 10,
                    light = FALSE)
  b <- a
  set.seed(310)
  full <- sim_fisheries(sim, a = a, b = b)$fishery
  a$lengths_cap <- 0
  a$reporting_cv <- 1
  a$sampling_coverage <- 0
  set.seed(310)
  reduced <- sim_fisheries(sim, a = a, b = b)$fishery
  expect_equal(full$joint_catch, reduced$joint_catch)
  expect_equal(full$mortality_ledger, reduced$mortality_ledger)
  columns <- c("operation", "fishery", "n", "cell", "effort")
  expect_equal(full$setdet[, ..columns], reduced$setdet[, ..columns])
})

test_that("existing custom mortality closures can return unnamed matrices", {
  withr::local_seed(311)
  sim <- sim_abundance(ages = 1:2, years = 2001:2002,
                        R = sim_R(log_mean = log(1000), log_sd = 0),
                        Z = function(ages, years) matrix(log(2), length(ages), length(years))) |>
    sim_distribution(grid = small_grid(), ays_covar = zero_covariance,
                       depth_par = function(data) rep(0, nrow(data)))
  expect_null(dimnames(sim$Z))
  spec <- sim_fishery(fishery_share = 0.1, catch_model = capture_all,
                       selectivity = function(age) 1, lengths_cap = 0)
  result <- sim_fisheries(sim, custom = spec)
  expect_equal(dimnames(result$Z), dimnames(result$N))
  expect_true(all(result$fishery$annual$status == "target_met"))
  expect_null(dimnames(sim$Z))
})
