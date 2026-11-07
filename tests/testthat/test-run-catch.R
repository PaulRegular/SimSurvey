catch_estimation_fixture <- function() {
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
  list(ages = 1:3, lengths = c(20, 30, 40), years = 2001, replicates = 1L,
       fishery = list(configurations = list(fleet = list()), setdet = operations, samp = sample))
}

catch_biology_fixture <- function(protocol = "measured") {
  sim <- catch_estimation_fixture()
  sample <- sim$fishery$samp
  sample$weight <- ifelse(sample$length == 20, 2, 8)
  sample$mature <- c(1, 0, 0, 1, 0, 0)
  sample$maturity_observed <- if (protocol == "measured") rep(TRUE, nrow(sample)) else sample$aged
  sample$mature[!sample$maturity_observed] <- NA_real_
  sample$maturity_probability <- if (protocol == "measured") rep(1, nrow(sample)) else
    sample$age_probability
  sim$fishery$samp <- sample
  sim
}

test_that("catch biological means use the stages needed for each observed trait", {
  measured <- run_catch(catch_biology_fixture())$fishery$catch_estimates
  aged <- run_catch(catch_biology_fixture("aged"))$fishery$catch_estimates
  expect_equal(measured$catch_at_length$mean_weight_kg, c(2, 8, NA))
  expect_equal(measured$catch_at_age$mean_weight_kg, c(2, 5, NA))
  expect_equal(measured$catch_at_length$prop_mature, c(2 / 7, 2 / 3, NA))
  expect_equal(aged$catch_at_length$prop_mature, c(3 / 4, 1, NA))
  expect_equal(measured$catch_at_age$prop_mature, c(1, 1 / 2, NA))
  expect_equal(aged$catch_at_age, measured$catch_at_age)
  expect_identical(aged$catch_at_length$weight_status, c("estimated", "estimated", "no_bin_sample"))
  expect_identical(aged$catch_at_age$maturity_status, c("estimated", "estimated", "no_bin_sample"))
  expect_equal(aged$catch_at_length$n_hat, c(140, 60, 0))
  expect_equal(aged$catch_at_age$n_hat, c(120, 80, 0))
})

test_that("catch biological support is separate from counts and other traits", {
  sim <- catch_biology_fixture()
  sample <- sim$fishery$samp
  sample$age_probability[] <- 0
  sample$aged[] <- FALSE
  sample$age[] <- NA_real_
  sim$fishery$samp <- sample
  result <- run_catch(sim)$fishery$catch_estimates
  expect_equal(result$catch_at_length$prop_mature, c(2 / 7, 2 / 3, NA))
  expect_true(all(result$catch_at_age$maturity_status == "no_aging_support"))

  sample$mature[] <- 0
  sample$weight[sample$length == 30] <- NA_real_
  sim$fishery$samp <- sample
  result <- run_catch(sim)$fishery$catch_estimates$catch_at_length
  expect_equal(result$prop_mature, c(0, 0, NA))
  expect_equal(result$mean_weight_kg, c(2, NA, NA))
  expect_identical(result$weight_status, c("estimated", "no_weight_support", "no_bin_sample"))
  expect_equal(result$n_hat, c(140, 60, 0))

  sample$maturity_probability[] <- 0
  sample$maturity_observed[] <- FALSE
  sample$mature[] <- NA_real_
  sim$fishery$samp <- sample
  result <- run_catch(sim)$fishery$catch_estimates$catch_at_length
  expect_true(all(is.na(result$prop_mature)))
  expect_identical(result$maturity_status, c("no_maturity_support", "no_maturity_support",
                                             "no_bin_sample"))
  expect_equal(result$n_hat, c(140, 60, 0))

  sim$fishery$setdet$reported_n[] <- 0
  result <- run_catch(sim)$fishery$catch_estimates$catch_at_length
  expect_true(all(is.na(result$mean_weight_kg) & is.na(result$prop_mature)))
  expect_true(all(result$weight_status == "zero_reported_catch" &
                    result$maturity_status == "zero_reported_catch" & result$n_hat == 0))
})

test_that("catch biology ignores latent values and retains length maturity when ages are missing", {
  sim <- catch_biology_fixture("aged")
  reference <- run_catch(sim)$fishery$catch_estimates
  sim$sim_weight <- function(x) stop("Latent weight model was used")
  sim$sim_maturity <- function(x) stop("Latent maturity model was used")
  sim$fishery$samp$true_age <- rep(3, 6)
  sim$fishery$samp$true_mature <- rep(1, 6)
  sim$fishery$samp$mature[!sim$fishery$samp$maturity_observed] <- 1
  expect_identical(run_catch(sim)$fishery$catch_estimates, reference)

  sim$fishery$samp$age[1] <- NA_real_
  result <- run_catch(sim)$fishery$catch_estimates
  expect_identical(result$catch_at_length, reference$catch_at_length)
  expect_true(all(result$catch_at_age$maturity_status == "no_aging_support"))

  sim <- catch_biology_fixture("aged")
  sim$fishery$samp$mature[1] <- NA_real_
  result <- run_catch(sim)$fishery$catch_estimates$catch_at_length
  expect_equal(result$prop_mature, c(NA, 1, NA))
  expect_identical(result$maturity_status, c("no_maturity_support", "estimated", "no_bin_sample"))
  expect_equal(result$n_hat, reference$catch_at_length$n_hat)

  sim <- catch_biology_fixture("aged")
  sim$fishery$samp[operation == 1 & length == 30, `:=`(aged = FALSE, age = NA_real_,
    age_probability = 1 / 3, maturity_probability = 1 / 3, maturity_observed = FALSE,
    mature = NA_real_)]
  result <- run_catch(sim)$fishery$catch_estimates$catch_at_length
  expect_identical(result$maturity_status, c("estimated", "no_bin_sample", "no_bin_sample"))
  expect_equal(result$prop_mature, c(3 / 4, NA, NA))
})

test_that("catch mean-weight truth divides aggregate kilograms by realized numbers", {
  sim <- catch_biology_fixture()
  sim$fishery$catch_at_age <- data.table(fishery = "fleet", sim = 1L, year = 2001,
                                        age = c(1, 2, 2, 3), n = c(50, 30, 20, 10),
                                        kg = c(100, 120, 140, 0))
  sim$fishery$catch_at_length <- data.table(fishery = "fleet", sim = 1L, year = 2001,
                                           length = c(20, 30, 40), n = c(100, 40, 10),
                                           kg = c(200, 280, NA))
  result <- sim |> run_catch() |> catch_error()
  expect_equal(result$fishery$catch_error$catch_at_age$mean_weight_true_kg, c(2, 5.2, 0))
  expect_equal(result$fishery$catch_error$catch_at_age$mean_weight_error_kg, c(0, -0.2, NA))
  expect_equal(result$fishery$catch_error$catch_at_length$mean_weight_true_kg, c(2, 7, NA))
  expect_equal(result$fishery$catch_error$catch_at_length$mean_weight_error_kg, c(0, 1, NA))
  stats <- result$fishery$catch_weight_error_stats
  expect_identical(key(stats$catch_at_length), "fishery")
  expect_equal(unlist(stats$catch_at_length[, -1]), c(ME = 0.5, MAE = 0.5, MSE = 0.5, RMSE = sqrt(0.5)))
  expect_equal(unlist(stats$catch_at_age[, -1]), c(ME = -0.1, MAE = 0.1, MSE = 0.02, RMSE = sqrt(0.02)))
  expect_null(run_catch(result)$fishery$catch_weight_error_stats)

  result$fishery$catch_at_age$kg <- NULL
  missing <- catch_error(result)$fishery
  expect_equal(missing$catch_error$catch_at_age$n_true, c(50, 50, 10))
  expect_true(all(is.na(missing$catch_error$catch_at_age$mean_weight_true_kg)))
  expect_true(all(is.na(unlist(missing$catch_weight_error_stats$catch_at_age[, -1]))))
  expect_identical(missing$catch_weight_error_stats$catch_at_length, stats$catch_at_length)
})

test_that("catch estimates use the corresponding inclusion stages and calibrate to reports", {
  sim <- catch_estimation_fixture()
  saved <- copy(sim)
  result <- run_catch(sim)$fishery
  estimates <- result$catch_estimates
  group <- c("fishery", "sim", "year")
  expect_identical(key(estimates$total_catch), group)
  expect_identical(key(estimates$catch_at_length), c(group, "length"))
  expect_identical(key(estimates$catch_at_age), c(group, "age"))
  expect_equal(estimates$total_catch$n_hat, 200)
  # Measured weights are 50 per fish in operation 1 and 25 in operation 2.
  # The short and long groups have weights 175 and 75 before calibration.
  expect_equal(estimates$catch_at_length$length, c(20, 30, 40))
  expect_equal(estimates$catch_at_length$n_hat, c(140, 60, 0))
  # Aged weights are 150, 50, and 50. Only the first is assigned to age 1.
  expect_equal(estimates$catch_at_age$age, 1:3)
  expect_equal(estimates$catch_at_age$n_hat, c(120, 80, 0))
  diagnostics <- result$catch_diagnostics
  expect_identical(diagnostics$length_status, "estimated")
  expect_identical(diagnostics$age_status, "estimated")
  expect_equal(diagnostics$reported_n, 200)
  expect_equal(diagnostics$n_operations, 3)
  expect_equal(diagnostics$n_positive_operations, 3)
  expect_equal(diagnostics$n_sampled_operations, 2)
  expect_equal(diagnostics$n_measured, 6)
  expect_equal(diagnostics$n_aged, 3)
  expect_equal(diagnostics$operation_coverage, 2 / 3)
  expect_equal(diagnostics$reported_catch_coverage, 140 / 200)
  expect_equal(sim, saved)
})

test_that("length estimates retain support when aging fails", {
  for (probability in c(0, NA_real_)) {
    sim <- catch_estimation_fixture()
    # A measured fish that could not be aged leaves the age composition unsupported.
    sim$fishery$samp$age_probability[2] <- probability
    result <- run_catch(sim)$fishery
    expect_identical(result$catch_diagnostics$length_status, "estimated")
    expect_identical(result$catch_diagnostics$age_status, "no_aging_support")
    expect_equal(result$catch_estimates$catch_at_length$n_hat, c(140, 60, 0))
    expect_true(all(is.na(result$catch_estimates$catch_at_age$n_hat)))
  }
  sim <- catch_estimation_fixture()
  sim$fishery$samp$aged[] <- FALSE
  sim$fishery$samp$age[] <- NA_real_
  result <- run_catch(sim)$fishery
  expect_identical(result$catch_diagnostics$length_status, "estimated")
  expect_identical(result$catch_diagnostics$age_status, "no_aging_support")
  expect_equal(result$catch_diagnostics$n_aged, 0)
  expect_equal(result$catch_estimates$catch_at_length$n_hat, c(140, 60, 0))
  expect_true(all(is.na(result$catch_estimates$catch_at_age$n_hat)))
})

test_that("operation and measurement support failures leave reported totals available", {
  for (probability in c(0, NA_real_)) {
    sim <- catch_estimation_fixture()
    sim$fishery$setdet$sampling_probability[3] <- probability
    result <- run_catch(sim)$fishery
    expect_identical(result$catch_diagnostics$length_status, "no_operation_support")
    expect_identical(result$catch_diagnostics$age_status, "no_operation_support")
    expect_equal(result$catch_estimates$total_catch$n_hat, 200)
    expect_true(all(is.na(result$catch_estimates$catch_at_length$n_hat)))
    expect_true(all(is.na(result$catch_estimates$catch_at_age$n_hat)))
  }
  sim <- catch_estimation_fixture()
  sim$fishery$samp <- sim$fishery$samp[operation == 2, ]
  result <- run_catch(sim)$fishery
  expect_identical(result$catch_diagnostics$length_status, "no_measurement_support")
  expect_identical(result$catch_diagnostics$age_status, "no_measurement_support")
  expect_equal(result$catch_estimates$total_catch$n_hat, 200)
  expect_true(all(is.na(result$catch_estimates$catch_at_length$n_hat)))
  expect_true(all(is.na(result$catch_estimates$catch_at_age$n_hat)))

  sim <- catch_estimation_fixture()
  sim$fishery$setdet$reported_n[1:2] <- 0
  result <- run_catch(sim)$fishery
  expect_equal(result$catch_diagnostics$n_measured, 6)
  expect_identical(result$catch_diagnostics$length_status, "no_measurement_support")
  expect_identical(result$catch_diagnostics$age_status, "no_measurement_support")
  expect_equal(result$catch_estimates$total_catch$n_hat, 60)
  expect_true(all(is.na(result$catch_estimates$catch_at_length$n_hat)))
  expect_true(all(is.na(result$catch_estimates$catch_at_age$n_hat)))
})

test_that("reported zero totals give zero catch estimates without sampling support", {
  sim <- catch_estimation_fixture()
  sim$fishery$setdet$n <- rep(50, 3)
  sim$fishery$setdet$reported_n[] <- 0
  sim$fishery$setdet$sampling_probability[] <- 0
  sim$fishery$samp$age_probability[] <- 0
  result <- catch_error(run_catch(sim))$fishery
  expect_identical(result$catch_diagnostics$length_status, "zero_reported_catch")
  expect_identical(result$catch_diagnostics$age_status, "zero_reported_catch")
  expect_true(is.na(result$catch_diagnostics$reported_catch_coverage))
  expect_equal(result$catch_estimates$total_catch$n_hat, 0)
  expect_equal(result$catch_estimates$catch_at_length$n_hat, rep(0, 3))
  expect_equal(result$catch_estimates$catch_at_age$n_hat, rep(0, 3))
  expect_equal(result$catch_error$total_catch$n_true, 150)
  expect_equal(result$catch_error$total_catch$error, -150)
})

test_that("repeated operation labels align by fleet, replicate, and year", {
  sim <- catch_estimation_fixture()
  operations <- sim$fishery$setdet
  sample <- sim$fishery$samp
  sim$replicates <- 1:2
  sim$years <- 2001:2002
  sim$fishery$configurations$other <- list()
  sim$fishery$setdet <- rbind(
    operations,
    copy(operations)[, `:=`(fishery = "other", sim = 2L, reported_n = reported_n * 2)],
    copy(operations)[, `:=`(year = 2002, reported_n = reported_n * 3)]
  )
  sim$fishery$samp <- rbind(sample,
                           copy(sample)[, `:=`(fishery = "other", sim = 2L)],
                           copy(sample)[, year := 2002])
  result <- run_catch(sim)$fishery
  estimates <- result$catch_estimates
  expect_equal(nrow(estimates$total_catch), 8)
  expect_equal(nrow(estimates$catch_at_length), 8 * 3)
  expect_equal(nrow(estimates$catch_at_age), 8 * 3)
  expect_equal(estimates$total_catch[fishery == "fleet" & sim == 1L, n_hat], c(200, 600))
  expect_equal(estimates$total_catch[fishery == "other" & sim == 2L, n_hat], c(400, 0))
  expect_equal(estimates$catch_at_length[fishery == "other" & sim == 2L & year == 2001, n_hat],
               c(280, 120, 0))
  expect_equal(estimates$catch_at_age[fishery == "fleet" & sim == 1L & year == 2002, n_hat],
               c(360, 240, 0))
  no_operations <- result$catch_diagnostics[n_operations == 0, ]
  expect_equal(nrow(no_operations), 5)
  expect_true(all(no_operations$length_status == "zero_reported_catch"))
  expect_true(all(no_operations$age_status == "zero_reported_catch"))
  expect_true(all(is.na(no_operations$operation_coverage)))
  expect_true(all(is.na(no_operations$reported_catch_coverage)))
  group <- c("fishery", "sim", "year")
  expect_true(all(estimates$catch_at_length[no_operations, on = group, n_hat] == 0))
  expect_true(all(estimates$catch_at_age[no_operations, on = group, n_hat] == 0))
})

test_that("catch estimation ignores latent fields and comparisons use independent truth tables", {
  sim <- catch_estimation_fixture()
  reference <- run_catch(sim)$fishery
  unavailable <- catch_error(run_catch(sim))$fishery$catch_error
  expect_true(all(vapply(unavailable, function(table) {
    all(table$truth_status == "unavailable" & is.na(table$n_true) & is.na(table$error))
  }, logical(1))))

  sim$fishery$annual <- data.table(fishery = "fleet", sim = 1L, year = 2001, reported_n = 600, n = 600)
  sim$fishery$setdet$n <- c(100, 40, 20)
  sim$fishery$setdet$n_measured <- rep(100L, 3)
  sim$fishery$samp$length_probability <- rep(0.01, 6)
  sim$fishery$samp$true_age <- rep(3, 6)
  sim$fishery$samp$operation_probability <- rep(0.01, 6)
  sim$fishery$samp$age[!sim$fishery$samp$aged] <- 99
  sim$fishery$samp <- rbind(sim$fishery$samp, copy(sim$fishery$samp[1, ])[
    , `:=`(measured = FALSE, aged = FALSE, length = 999, age = 99, age_probability = NA_real_)])
  sim$fishery$catch_at_age <- data.table(fishery = "fleet", sim = 1L, year = 2001,
                                        operation = c(1, 1, 2, 3), age = c(1, 2, 2, 2),
                                        n = c(80, 20, 40, 20))
  sim$fishery$catch_at_length <- data.table(fishery = "fleet", sim = 1L, year = 2001,
                                           operation = c(1, 1, 2, 2, 3),
                                           length = c(20, 30, 20, 30, 30),
                                           n = c(70, 30, 20, 20, 20))
  result <- run_catch(sim) |> catch_error()
  expect_identical(result$fishery$catch_estimates, reference$catch_estimates)
  expect_identical(result$fishery$catch_diagnostics, reference$catch_diagnostics)
  comparisons <- result$fishery$catch_error
  expect_equal(comparisons$total_catch$n_true, 160)
  expect_equal(comparisons$total_catch$error, 40)
  expect_equal(comparisons$catch_at_length$n_true, c(90, 70, 0))
  expect_equal(comparisons$catch_at_length$error, c(50, -10, 0))
  expect_equal(comparisons$catch_at_age$n_true, c(80, 80, 0))
  expect_equal(comparisons$catch_at_age$error, c(40, 0, 0))
  for (field in names(comparisons)) {
    missing_truth <- copy(result)
    if (field == "total_catch") missing_truth$fishery$setdet$n <- NULL else
      missing_truth$fishery[[field]] <- NULL
    compared <- catch_error(missing_truth)$fishery$catch_error
    expect_true(all(is.na(compared[[field]]$n_true)))
    expect_true(all(is.na(compared[[field]]$error)))
    expect_true(all(compared[[field]]$truth_status == "unavailable"))
    expect_identical(compared[setdiff(names(compared), field)],
                     comparisons[setdiff(names(comparisons), field)])
  }
  expect_null(run_catch(result)$fishery$catch_error)
  expect_error(catch_error(sim), "Run run_catch")
})

test_that("a monitored census matches realized catch in all tables", {
  withr::local_seed(224)
  config <- sim_fishery(target = 0.25, catch_model = capture_all, selectivity = function(x) 1,
                        subset_cells = cell == 1, sampling_coverage = 1,
                        lengths_cap = 1000, ages_cap = 1000)
  empty <- config
  empty$target <- 0
  pop <- fishery_population(years = 2001:2002)
  pop$lengths <- as.numeric(rownames(pop$N_at_length))
  sim <- pop |>
    sim_fisheries(a = config, b = config, empty = empty, n_sims = 2) |>
    run_catch() |>
    catch_error()
  expect_equal(nrow(sim$fishery$catch_estimates$total_catch), 3 * 2 * 2)
  expect_equal(nrow(sim$fishery$catch_estimates$catch_at_length), 3 * 2 * 2 * 2)
  expect_equal(nrow(sim$fishery$catch_estimates$catch_at_age), 3 * 2 * 2 * 2)
  for (comparison in sim$fishery$catch_error) {
    expect_equal(comparison$n_hat, comparison$n_true)
    expect_equal(comparison$error, rep(0, nrow(comparison)))
    expect_true(all(comparison$truth_status == "available"))
  }
  for (trait in c("catch_at_age", "catch_at_length")) {
    comparison <- sim$fishery$catch_error[[trait]]
    expect_equal(comparison$mean_weight_kg, comparison$mean_weight_true_kg)
    errors <- comparison$mean_weight_error_kg[!is.na(comparison$mean_weight_error_kg)]
    expect_equal(errors, numeric(length(errors)))
  }
  expect_true(all(sim$fishery$catch_diagnostics[fishery == "empty", length_status] ==
                    "zero_reported_catch"))
  expect_true(all(sim$fishery$catch_diagnostics[fishery == "empty", age_status] ==
                    "zero_reported_catch"))
  no_operations <- pop |> sim_fisheries(empty = empty) |> run_catch() |> catch_error()
  for (comparison in no_operations$fishery$catch_error) {
    expect_true(all(comparison$n_hat == 0 & comparison$n_true == 0))
    expect_true(all(comparison$truth_status == "available"))
  }
})

test_that("catch estimation rejects malformed reports and biological sample keys", {
  sim <- catch_estimation_fixture()
  sim$fishery$setdet$reported_n[1] <- NA_real_
  expect_error(run_catch(sim), "finite nonnegative")
  sim <- catch_estimation_fixture()
  sim$fishery$setdet <- rbind(sim$fishery$setdet, sim$fishery$setdet[1, ])
  expect_error(run_catch(sim), "operation keys must be unique")
  sim <- catch_estimation_fixture()
  sim$fishery$samp$operation[1] <- 99L
  expect_error(run_catch(sim), "must match selected operations")
  sim <- catch_estimation_fixture()
  sim$fishery$samp$length[1] <- 99
  expect_error(run_catch(sim), "native simulation length groups")
  sim <- catch_estimation_fixture()
  sim$fishery$samp$measured[1] <- FALSE
  expect_error(run_catch(sim), "Aged fish must be measured")
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
