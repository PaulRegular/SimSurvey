test_that("shortfalls retain analytical reference and realized mortality", {
  withr::local_seed(711)
  population <- fishery_population(100)
  population$sp_N$N <- ifelse(population$sp_N$cell == 1, 5,
                             ifelse(population$sp_N$cell == 2, 95, 0))
  result <- sim_fishery(population, 0.4, capture_all, selectivity = function(age) 1,
                         subset_cells = cell == 1)$fishery
  reference <- result$reference
  expect_equal(reference$N, 100)
  expect_equal(reference$Z, log(2))
  expect_equal(reference$death_budget, 50)
  expect_equal(reference$total_deaths, 50)
  expect_identical(reference$total_deaths, reference$death_budget)
  expect_equal(reference$catch_reference, 20)
  expect_equal(reference$F_reference, 0.4 * log(2))
  expect_equal(reference$catch_realized, 5)
  expect_equal(reference$U_realized, 0.05)
  expect_equal(reference$F_realized, 0.1 * log(2))
  expect_equal(result$annual$shortfall, 15)
  expect_equal(result$annual$status, "exhausted")
  expect_equal(sum(result$mortality_ledger$other_deaths), 45)
  expect_equal(reference$total_deaths, sum(result$mortality_ledger$fishery_deaths +
                                            result$mortality_ledger$other_deaths))
  expect_equal(sum(result$mortality_ledger$survivors), 50)
})

test_that("mortality diagnostics retain zeros, undefined rates and rounding", {
  withr::local_seed(712)
  population <- fishery_population(c(100, 0))
  result <- fishery_run(population)$fishery$reference
  expect_equal(result$catch_realized, c(25, 0))
  expect_equal(result$U_realized, c(0.25, NA))
  expect_equal(result$F_realized, c(log(2) / 2, NA))
  population$Z[] <- 0
  result <- fishery_run(population)$fishery$reference
  expect_equal(result$catch_realized, c(0, 0))
  expect_equal(result$F_realized, c(0, NA))
  zero <- sim_fishery(fishery_population(), 0, capture_all)$fishery$reference
  expect_equal(zero$F_realized, c(0, 0))
  rounded <- sim_fishery(fishery_population(3), 1, capture_all,
                          selectivity = function(age) 1)$fishery$reference
  expect_equal(rounded$catch_reference, 1.5)
  expect_equal(rounded$catch_realized, 2)
  expect_equal(rounded$total_deaths, 2)
  expect_equal(rounded$F_realized, 4 * log(2) / 3)
  expect_gt(rounded$F_realized, rounded$Z)
})

test_that("realized F follows age composition and true catches before reporting", {
  withr::local_seed(713)
  population <- fishery_population()
  capture_young <- function(available, operation) available$N * (available$age == 1)
  set.seed(714)
  full <- sim_fishery(population, 0.4, capture_young, selectivity = function(age) 1,
                       reporting_factor = 1)$fishery
  set.seed(714)
  unreported <- sim_fishery(population, 0.4, capture_young, selectivity = function(age) 1,
                             reporting_factor = 0)$fishery
  expect_equal(full$reference, unreported$reference)
  expect_equal(full$annual$n, full$annual$target)
  expect_equal(full$reference$catch_realized, c(40, 0))
  expect_equal(full$reference$F_realized, c(0.8 * log(2), 0))
  expect_gt(full$reference$F_realized[1], full$reference$F_reference[1])
  expect_equal(sum(full$reference$catch_realized), sum(full$annual$n))
})
