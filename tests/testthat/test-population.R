test_that("recruitment closures implement constant and random-walk recruitment", {
  withr::local_seed(438)
  years <- 2001:2003
  expect_equal(unname(sim_R(log_mean = log(c(100, 200, 300)), log_sd = 0)(years)),
               c(100, 200, 300))
  expect_equal(names(sim_R(log_sd = 0)(years)), as.character(years))

  set.seed(438)
  expected <- exp(log(100) + cumsum(stats::rnorm(3, sd = 0.2)))
  set.seed(438)
  expect_equal(unname(sim_R(log_mean = log(100), log_sd = 0.2)(years)), expected)
  expect_error(sim_R(log_mean = c(1, 2))(years), "number of log_means")
})

test_that("mortality and starting abundance retain age/year contracts", {
  withr::local_seed(438)
  mortality <- sim_Z(log_mean = log(0.2), log_sd = 0)(ages = 2:4, years = 2001:2003)
  expect_equal(mortality, matrix(0.2, 3, 3,
                                dimnames = list(age = 2:4, year = 2001:2003)))
  expect_equal(sim_N0()(R0 = 100, Z0 = rep(0.2, 3), ages = 2:4),
               100 * exp(-0.2 * 0:2))
  expect_equal(sim_N0(c(80, 50))(R0 = 100, Z0 = rep(0.2, 3), ages = 2:4),
               c(100, 80, 50))
  expect_error(sim_Z(log_mean = matrix(0, 2, 2))(ages = 1:3, years = 1:3),
               "matrix of log means")
})

test_that("cohort survival agrees with an independently specified population", {
  withr::local_seed(438)
  sim <- sim_abundance(
    ages = 2:4, years = 2001:2003,
    R = sim_R(log_mean = log(c(100, 200, 300)), log_sd = 0),
    Z = sim_Z(log_mean = log(0.2), log_sd = 0),
    N0 = sim_N0(c(80, 50))
  )
  expected <- rbind(c(100, 200, 300),
                    c(80, 100 * exp(-0.2), 200 * exp(-0.2)),
                    c(50, 80 * exp(-0.2), 100 * exp(-0.4)))
  dimnames(expected) <- list(age = 2:4, year = 2001:2003)
  expect_equal(sim$N, expected)
  expect_equal(colSums(sim$N_at_length), colSums(expected))
  expect_equal(sim$N0, expected[, 1], ignore_attr = TRUE)
  expect_true(is.function(sim$sim_length))
})

test_that("growth probabilities and DFO length grouping are consistent", {
  withr::local_seed(438)
  growth <- sim_vonB(Linf = 80, L0 = 10, K = 0.2, log_sd = 0.1, length_group = 3)
  key <- growth(age = 1:3, length_age_key = TRUE)
  expect_true(is.matrix(key))
  expect_equal(unname(colSums(key)), rep(1, 3), tolerance = 1e-8)
  expect_true(all(key >= 0 & key <= 1))
  exact_growth <- sim_vonB(Linf = 80, L0 = 10, K = 0.2, log_sd = 0, digits = 3)
  expect_equal(exact_growth(1:3), round(80 - 70 * exp(-0.2 * 1:3), 3))
  expect_equal(group_lengths(c(0, 1, 2, 3, 5, 6), 3), c(1, 1, 1, 4, 4, 7))
  expect_equal(group_lengths(c(0, 1, 2, 3), 1), 0:3)
  expect_equal(group_lengths(c(0, 0.5, 1, 1.5), 0.5), c(0, 0.5, 1, 1.5))
})

test_that("age-to-length conversion matches matrix multiplication", {
  abundance <- matrix(c(100, 40, 80, 20), 2,
                       dimnames = list(age = 1:2, year = 2001:2002))
  key <- matrix(c(0.8, 0.2, 0.25, 0.75), 2,
                dimnames = list(length = c(10, 20), age = 1:2))
  result <- convert_N(abundance, key)
  expect_equal(result, key %*% abundance)
  expect_equal(colSums(result), colSums(abundance))
})
