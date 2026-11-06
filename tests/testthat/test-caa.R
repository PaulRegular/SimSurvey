caa_fixture <- function() {
  operations <- data.table(fishery = "fleet", sim = 1L, year = 2001,
                           operation = 1:3, reported_n = c(100, 40, 60),
                           sampling_probability = c(0.5, 0.8, 0.25),
                           sampled = c(TRUE, TRUE, FALSE))
  # Operation 1 ages one of three short fish and its one long fish.
  # Operation 2 ages one of two fish at random.
  sample <- data.table(fishery = "fleet", sim = 1L, year = 2001,
                       operation = rep(1:2, c(4, 2)), measured = TRUE,
                       length = c(20, 20, 20, 30, 20, 30),
                       aged = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE),
                       age = c(1, NA, NA, 2, 2, NA),
                       age_probability = c(rep(1 / 3, 3), 1, 0.5, 0.5))
  list(ages = 1:3, years = 2001, replicates = 1L,
       fishery = list(configurations = list(fleet = list()), setdet = operations, samp = sample))
}

test_that("catch-at-age weights cover both sampling stages and calibrate to all reports", {
  sim <- caa_fixture()
  result <- run_caa(sim)$fishery
  # The aged weights are 150, 50, and 50. Only the first is assigned to age 1.
  expect_equal(result$caa$age, 1:3)
  expect_equal(result$caa$n_hat, c(120, 80, 0))
  expect_equal(result$caa_diagnostics$status, "estimated")
  expect_equal(result$caa_diagnostics$reported_n, 200)
  expect_equal(result$caa_diagnostics$n_operations, 3)
  expect_equal(result$caa_diagnostics$n_positive_operations, 3)
  expect_equal(result$caa_diagnostics$n_sampled_operations, 2)
  expect_equal(result$caa_diagnostics$n_measured, 6)
  expect_equal(result$caa_diagnostics$n_aged, 3)
  expect_equal(result$caa_diagnostics$operation_coverage, 2 / 3)
  expect_equal(result$caa_diagnostics$reported_catch_coverage, 140 / 200)
  expect_null(sim$fishery$caa)
  expect_false("n_measured" %in% names(sim$fishery$setdet))
})

test_that("catch-at-age distinguishes absent support from reported zeros", {
  for (pi in c(0, NA_real_)) {
    sim <- caa_fixture()
    sim$fishery$setdet$sampling_probability[3] <- pi
    result <- run_caa(sim)$fishery
    expect_equal(result$caa_diagnostics$status, "no_operation_support")
    expect_true(all(is.na(result$caa$n_hat)))
  }

  sim <- caa_fixture()
  sim$fishery$samp <- sim$fishery$samp[operation == 2, ]
  result <- run_caa(sim)$fishery
  expect_equal(result$caa_diagnostics$status, "no_measurement_support")
  expect_true(all(is.na(result$caa$n_hat)))

  sim <- caa_fixture()
  sim$fishery$samp$age_probability[2] <- 0
  result <- run_caa(sim)$fishery
  expect_equal(result$caa_diagnostics$status, "no_aging_support")
  expect_true(all(is.na(result$caa$n_hat)))

  sim$fishery$setdet$sampled[] <- FALSE
  sim$fishery$samp <- sim$fishery$samp[0, ]
  result <- run_caa(sim)$fishery
  expect_equal(result$caa_diagnostics$status, "no_aging_support")
  expect_true(all(is.na(result$caa$n_hat)))

  sim <- caa_fixture()
  sim$fishery$setdet$n <- rep(50, 3)
  sim$fishery$setdet$reported_n[] <- 0
  sim$fishery$setdet$sampling_probability[] <- 0
  sim$fishery$samp$age_probability[] <- 0
  result <- run_caa(sim)$fishery
  expect_equal(result$caa_diagnostics$status, "zero_reported_catch")
  expect_equal(result$caa$n_hat, rep(0, 3))
  expect_true(is.na(result$caa_diagnostics$reported_catch_coverage))
})

test_that("catch-at-age aligns repeated operation labels across fleets and replicates", {
  sim <- caa_fixture()
  sim$replicates <- 1:2
  sim$fishery$configurations$other <- list()
  sim$fishery$setdet <- rbind(sim$fishery$setdet, copy(sim$fishery$setdet))
  sim$fishery$setdet[4:6, `:=`(fishery = "other", sim = 2L, reported_n = reported_n * 2)]
  sim$fishery$samp <- rbind(sim$fishery$samp, copy(sim$fishery$samp))
  sim$fishery$samp[7:12, `:=`(fishery = "other", sim = 2L)]
  result <- run_caa(sim)$fishery
  expect_equal(result$caa[fishery == "fleet" & sim == 1L, n_hat], c(120, 80, 0))
  expect_equal(result$caa[fishery == "other" & sim == 2L, n_hat], c(240, 160, 0))
  expect_equal(result$caa[fishery == "fleet" & sim == 2L, n_hat], rep(0, 3))
  expect_equal(result$caa[fishery == "other" & sim == 1L, n_hat], rep(0, 3))
  no_operations <- result$caa_diagnostics[result$caa_diagnostics$n_operations == 0, ]
  expect_equal(no_operations$status, rep("zero_reported_catch", 2))
  expect_true(all(is.na(no_operations$operation_coverage)))
  expect_true(all(is.na(no_operations$reported_catch_coverage)))
})

test_that("catch-at-age uses no hidden fields and the comparison joins truth separately", {
  sim <- caa_fixture()
  reference <- run_caa(sim)$fishery
  unavailable <- caa_error(run_caa(sim))$fishery$caa_error
  expect_equal(unavailable$truth_status, rep("unavailable", 3))
  expect_true(all(is.na(unavailable$n_true)))
  expect_true(all(is.na(unavailable$error)))

  sim$fishery$annual <- data.table(fishery = "fleet", sim = 1L, year = 2001, reported_n = 600, n = 600)
  sim$fishery$setdet$n <- rep(200, 3)
  sim$fishery$samp$length_probability <- rep(0.01, 6)
  sim$fishery$samp$true_age <- rep(3, 6)
  sim$fishery$samp$operation_probability <- rep(0.01, 6)
  sim$fishery$samp$age[!sim$fishery$samp$aged] <- 3
  sim$fishery$catch_at_age <- data.table(fishery = "fleet", sim = 1L, year = 2001,
                                        operation = c(1, 1, 2), age = c(1, 2, 2), n = c(80, 20, 40))
  result <- run_caa(sim) |> caa_error()
  expect_identical(result$fishery$caa, reference$caa)
  expect_identical(result$fishery$caa_diagnostics, reference$caa_diagnostics)
  expect_equal(result$fishery$caa_error$n_true, c(80, 60, 0))
  expect_equal(result$fishery$caa_error$error, c(40, 20, 0))
  expect_equal(result$fishery$caa_error$truth_status, rep("available", 3))
  expect_null(run_caa(result)$fishery$caa_error)
  expect_error(caa_error(sim), "Run run_caa")
})

test_that("a monitored census matches realized catch and empty fisheries retain zeros", {
  withr::local_seed(224)
  config <- sim_fishery(target = 0.25, catch_model = capture_all, selectivity = function(x) 1,
                        subset_cells = cell == 1, sampling_coverage = 1,
                        lengths_cap = 1000, ages_cap = 1000)
  empty <- config
  empty$target <- 0
  sim <- fishery_population(years = 2001:2002) |>
    sim_fisheries(a = config, b = config, empty = empty, n_sims = 2) |>
    run_caa() |>
    caa_error()
  expect_equal(nrow(sim$fishery$caa), 3 * 2 * 2 * 2)
  expect_equal(sim$fishery$caa_error$n_hat, sim$fishery$caa_error$n_true)
  expect_equal(sim$fishery$caa_error$error, rep(0, nrow(sim$fishery$caa)))
  expect_equal(sim$fishery$caa_diagnostics[fishery == "empty", status], rep("zero_reported_catch", 4))

  no_operations <- fishery_population() |>
    sim_fisheries(empty = empty) |>
    run_caa() |>
    caa_error()
  expect_equal(no_operations$fishery$caa_error$n_true, c(0, 0))
  expect_equal(no_operations$fishery$caa_error$n_hat, c(0, 0))
  expect_equal(no_operations$fishery$caa_diagnostics$status, "zero_reported_catch")
})

test_that("catch-at-age rejects malformed operation reports and sampling keys", {
  sim <- caa_fixture()
  sim$fishery$setdet$reported_n[1] <- NA_real_
  expect_error(run_caa(sim), "finite nonnegative")
  sim <- caa_fixture()
  sim$fishery$setdet <- rbind(sim$fishery$setdet, sim$fishery$setdet[1, ])
  expect_error(run_caa(sim), "operation keys must be unique")
  sim <- caa_fixture()
  sim$fishery$samp$operation[1] <- 99L
  expect_error(run_caa(sim), "must match selected operations")
})

test_that("grouped monitoring counts preserve operation order, numeric types and zeros", {
  sample <- data.table(operation = c(2L, 1L, 2L, 1L, 2L), aged = c(TRUE, TRUE, FALSE, FALSE, TRUE),
                       maturity_observed = c(TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_identical(.fishery_sample_counts(c(3L, 1L, 2L), sample),
                   list(n_measured = c(0, 2, 3), n_aged = c(0, 1, 2), n_maturity = c(0, 1, 1)))
  expect_identical(.fishery_sample_counts(c(1L, 2L), sample[0, ]),
                   list(n_measured = c(0, 0), n_aged = c(0, 0), n_maturity = c(0, 0)))
  expect_identical(.fishery_sample_counts(integer(), sample[0, ]),
                   list(n_measured = numeric(), n_aged = numeric(), n_maturity = numeric()))
})
