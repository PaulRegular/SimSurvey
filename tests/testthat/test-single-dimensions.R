test_that("single-age length keys retain their matrix dimensions", {
  withr::local_seed(438)
  key <- sim_vonB()(age = 1, length_age_key = TRUE)
  expect_true(is.matrix(key))
  expect_equal(ncol(key), 1L)
  expect_equal(unname(colSums(key)), 1, tolerance = 1e-8)
  sim <- sim_abundance(ages = 1, years = 2001:2002)
  expect_equal(dim(sim$N), c(1L, 2L))
})

test_that("rounding retains one-age and one-year population matrices", {
  withr::local_seed(438)
  for (ages in list(1, 1:3)) {
    for (years in list(2001, 2001:2002)) {
      sim <- sim_abundance(ages = ages, years = years) |>
        sim_distribution(grid = small_grid(), ays_covar = zero_covariance,
                          depth_par = function(data) rep(0, nrow(data))) |>
        round_sim()
      expect_equal(dim(sim$N), c(length(ages), length(years)))
      expect_equal(colSums(sim$N), colSums(sim$N_at_length))
      expect_equal(sum(sim$N), sum(sim$sp_N$N))
    }
  }
})
