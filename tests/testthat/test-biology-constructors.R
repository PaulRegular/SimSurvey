test_that("biological constructors preserve units, aliases and RNG", {
  withr::local_seed(401)
  rng <- .Random.seed
  weight <- sim_weight(alpha = 0.001, beta = 3)
  maturity <- sim_maturity()
  growth <- sim_growth()
  expect_identical(.Random.seed, rng)
  expect_identical(sim_growth, sim_vonB)
  expect_equal(weight(c(0, 10, 20)), c(0, 1, 8))
  expect_equal(sim_weight()(1), 10^-5.2106)
  expect_equal(maturity(47.81), 0.5)
  expect_equal(maturity(c(30, 47.81, 60)), sim_logistic(k = 0.28, x50 = 47.81)(c(30, 47.81, 60)))
  expect_equal(sim_logistic()(1:5), 1 / (1 + exp(-2 * ((1:5) - 3))))
  expect_equal(sim_logistic(2, 3)(c(-1e6, 1e6)), c(0, 1))
  expect_error(sim_weight(-1), "alpha")
  expect_error(sim_weight(beta = 0), "beta")
  expect_error(weight(-1), "Lengths")
  expect_error(sim_maturity(k = -1), "steepness")
  expect_error(sim_logistic(x50 = Inf), "finite")
})

test_that("x0 is a deprecated named alias, without changing positional calls", {
  withr::local_options(lifecycle_verbosity = "warning")
  expect_warning(old <- sim_logistic(k = 0.5, x0 = 4), class = "lifecycle_warning_deprecated")
  expect_equal(old(1:8), sim_logistic(0.5, 4)(1:8))
  expect_error(sim_logistic(x50 = 4, x0 = 3), "both")
})

