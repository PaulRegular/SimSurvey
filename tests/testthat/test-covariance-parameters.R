test_that("age-specific innovation SD scales independent age fields", {
  withr::local_seed(438)
  x <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  simulate <- function(sd) {
    sim_ays_covar(sd = sd, phi_age = 0, phi_year = 0, group_ages = NULL)(
      x, ages = 1:3, years = 1:3, cells = 1:3
    )
  }
  baseline <- simulate(1)
  set.seed(438)
  result <- simulate(c(1, 2, 5))
  expect_equal(result, sweep(baseline, 1, c(1, 2, 5), "*"), tolerance = 1e-12)
})

test_that("varying AR parameters use the current innovation scale", {
  withr::local_seed(438)
  x <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)
  sd <- c(0.1, 0.2, 0.5)
  pa <- c(0.1, 0.3, 0.7)
  py <- c(0.2, 0.4, 0.6)
  noise <- array(stats::rnorm(18), dim = c(2, 3, 3))
  spatial <- t(chol(exp(-as.matrix(stats::dist(x)) / 10)))
  expected <- array(0, c(3, 3, 2),
                     dimnames = list(age = as.character(1:3),
                                     year = as.character(1:3), cell = as.character(1:2)))
  for (y in 1:3) {
    for (a in 1:3) {
      mean <- if (y > 1) py[y] * expected[a, y - 1, ] else 0
      if (a > 1) {
        previous <- if (y > 1) py[y] * expected[a - 1, y - 1, ] else 0
        mean <- mean + pa[a] * (expected[a - 1, y, ] - previous)
      }
      scale <- sd[a] / (if (a == 1) sqrt(1 - pa[a] ^ 2) else 1) /
        (if (y == 1) sqrt(1 - py[y] ^ 2) else 1)
      expected[a, y, ] <- mean + scale * (spatial %*% noise[, a, y])
    }
  }
  set.seed(438)
  result <- sim_ays_covar(sd = sd, phi_age = pa, phi_year = py, group_ages = NULL,
                          range = 10, model = "exponential")(x, 1:3, 1:3, 1:2)
  expect_equal(result, expected, tolerance = 1e-12)
})
