test_that("reference targets, mortality and retained effort balance analytically", {
  withr::local_seed(21)
  sim <- fishery_population()
  saved <- data.table::copy(sim)
  result <- fishery_run(sim, gear_amount = 3, duration = 4, weight = function(length) length / 10,
                         lengths_cap = 8, ages_cap = 2, light = FALSE)
  f <- result$fishery
  expect_equal(f$reference$F_reference, rep(log(2) / 2, 2))
  expect_equal(f$reference$catch_reference, c(25, 25))
  expect_equal(f$annual$target, 50)
  expect_equal(f$annual$n, 50)
  expect_equal(f$annual$status, "target_met")
  expect_equal(f$setdet$effort, 12)
  expect_true(f$setdet$quota_limited)
  expect_true(f$setdet$mortality_limited)
  expect_equal(f$annual$cpue, 50 / 12)
  expect_equal(sum(f$joint_catch$n), 50)
  expect_equal(f$setdet$kg, sum(f$joint_catch$n * f$joint_catch$length / 10))
  ledger <- f$mortality_ledger
  expect_equal(ledger$starting, ledger$fishery_deaths + ledger$other_deaths + ledger$survivors)
  deaths <- ledger[, list(deaths = sum(fishery_deaths + other_deaths)), by = "age"]
  expect_equal(deaths$deaths, c(50, 50))
  expect_true(all(ledger$survivors >= 0))
  expect_equal(nrow(f$samp), 8)
  expect_true(all(is.na(f$samp$age[!f$samp$aged])))
  expect_true(all(f$samp$length_probability == 8 / 50))
  expect_equal(sim, saved)
  expect_equal(result$N, sim$N)
})

test_that("empty and unreachable fisheries return explicit statuses", {
  withr::local_seed(22)
  sim <- fishery_population()
  zero <- sim_fishery(sim, 0, capture_all, selectivity = function(x) 1)
  expect_equal(zero$fishery$annual$target, 0)
  expect_equal(nrow(zero$fishery$setdet), 0)
  expect_true(is.na(zero$fishery$annual$cpue))
  sim$Z[] <- 0
  expect_equal(fishery_run(sim)$fishery$annual$n, 0)
  sim <- fishery_population()
  inaccessible <- sim_fishery(sim, 0.5, capture_all, selectivity = function(x) 1,
                              subset_cells = FALSE)
  expect_equal(inaccessible$fishery$annual$status, "exhausted")
  expect_equal(inaccessible$fishery$annual$shortfall, 50)
  unsuccessful <- sim_fishery(sim, 0.5, function(available, operation) rep(0, nrow(available)),
                               selectivity = function(x) 1, max_operations = 3)
  expect_equal(unsuccessful$fishery$setdet$n, rep(0, 3))
  expect_equal(unsuccessful$fishery$setdet$effort, rep(1, 3))
  expect_equal(unsuccessful$fishery$annual$status, "operation_limit")
  expect_equal(unsuccessful$fishery$annual$cpue, 0)
  no_selection <- sim_fishery(sim, 1, capture_all, selectivity = function(x) 0)
  expect_equal(no_selection$fishery$annual$target, 0)
})

test_that("capture is bounded, responds to area and obeys capacity", {
  withr::local_seed(23)
  available <- data.frame(age = 1:3, length = 10:12, N = c(100, 200, 300),
                          selectivity = c(1, 0, 1), cell_area = 1)
  op <- list(gear_amount = 3, duration = 2, effort = 6)
  expect_equal(sim_capture(0)(available, op), c(0, 0, 0))
  expect_equal(sim_capture(1e8)(available, op), c(100, 0, 300))
  caught <- sim_capture(1e8, capacity_per_gear = 2.5)(available, op)
  expect_equal(sum(caught), 7)
  expect_equal(caught[2], 0)
  expect_equal(sim_capture(1, 0)(available, op), c(0, 0, 0))
  expect_error(sim_capture(-1), "q must")
  expect_error(sim_capture(1, NA), "capacity")
  # q = log(2), G*T = 6, A = 6 gives exactly one-half encounter probability.
  available$cell_area <- 6
  set.seed(25)
  expected <- stats::rbinom(3, c(100, 200, 300), c(0.5, 0, 0.5))
  set.seed(25)
  expect_equal(sim_capture(log(2))(available, op), expected)
})

test_that("callback and control domains are enforced", {
  withr::local_seed(24)
  sim <- fishery_population()
  expect_error(sim_fishery(sim, 1.1, capture_all), "fishery_share")
  expect_error(sim_fishery(sim, 0.2, capture_all, selectivity = function(x) 2), "probabilities")
  expect_error(fishery_run(sim, gear_amount = function(operation) -1), "gear_amount")
  expect_error(fishery_run(sim, duration = 0), "duration")
  expect_error(sim_fishery(sim, 0.2, function(available, operation) available$N + 1), "whole-fish")
  expect_error(sim_fishery(sim, 0.2, function(available, operation) 0.5), "whole-fish")
  expect_error(sim_fishery(sim, 0.2, capture_all, subset_cells = NA), "subset_cells")
})

test_that("configuration creation consumes no random numbers", {
  withr::local_seed(26)
  before <- .Random.seed
  spec <- sim_fishery(fishery_share = 0.3, catch_model = sim_capture(0.1),
                       subset_cells = depth > 20)
  expect_s3_class(spec, "sim_fishery_spec")
  expect_identical(.Random.seed, before)
  expect_identical(spec$subset_cells, quote(depth > 20))
})

test_that("large stocks stay aggregate and singleton experiments work", {
  withr::local_seed(27)
  result <- fishery_run(fishery_population(1e10), max_operations = 2,
                         sampling_coverage = 0, light = FALSE)
  expect_equal(result$fishery$annual$n, 2.5e9)
  expect_equal(nrow(result$fishery$joint_catch), 1)
  expect_equal(nrow(result$fishery$mortality_ledger), 16)
  expect_equal(nrow(result$fishery$samp), 0)
})

test_that("integer death budgets preserve the rounded total", {
  withr::local_seed(28)
  sim <- fishery_population(c(3, 3))
  result <- sim_fishery(sim, 0, capture_all)
  ledger <- result$fishery$mortality_ledger
  expect_equal(sum(ledger$other_deaths), 3)
  expect_equal(sum(ledger$survivors), 3)
})

test_that("repeated operations deplete stocks and retain unsuccessful effort", {
  withr::local_seed(29)
  seen <- numeric()
  capture <- function(available, operation) {
    seen <<- c(seen, sum(available$N))
    .thin_counts(available$N, 7)
  }
  result <- sim_fishery(fishery_population(100), 0.5, capture,
                         selectivity = function(x) 1, subset_cells = cell == 1,
                         sampling_coverage = 0)$fishery
  expect_equal(result$setdet$n, c(7, 7, 7, 4))
  expect_equal(seen, c(100, 93, 86, 79))
  expect_equal(result$setdet$quota_limited, c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(result$annual$effort, 4)
  # Restrict access to a small portion of the initial population.
  sim <- fishery_population(100)
  sim$sp_N$N[sim$sp_N$cell == 1] <- 10
  sim$sp_N$N[sim$sp_N$cell == 2] <- 90
  result <- sim_fishery(sim, 1, capture_all, selectivity = function(x) 1,
                         subset_cells = cell == 1)$fishery
  expect_equal(result$annual$n, 10)
  expect_equal(result$annual$target, 50)
  expect_equal(result$annual$status, "exhausted")
})

test_that("location probabilities separate density, knowledge and exploration", {
  expect_equal(.location_weights(c(1, 3), c(0, 0), 2, 0.2), c(0.18, 0.82))
  expect_equal(.location_weights(c(0, 10), c(-2, 100), 0, 0), c(0.5, 0.5))
  expect_equal(.location_weights(c(1, 3), c(-2, -4), 2, 0), c(0.5, 0.5))
  expect_equal(.location_weights(c(1, 3), c(2, -2), 1, 0), c(0.75, 0.25))
})
