test_that("the small grid has unique cells, expected resolution, and one stratum", {
  grid <- small_grid()
  cells <- as.data.frame(grid)
  expect_s3_class(grid, "stars")
  expect_equal(nrow(cells), 16)
  expect_equal(anyDuplicated(cells$cell), 0L)
  expect_equal(length(unique(cells$strat)), 1L)
  expect_equal(abs(stars::st_res(grid)), c(x = 5, y = 5))
})

test_that("distribution conserves abundance and aligns age/year/cell indices", {
  withr::local_seed(438)
  sim <- small_population()
  totals <- with(sim$sp_N, tapply(N, list(age, year), sum))
  expect_equal(unname(totals), unname(sim$N))
  expect_equal(sort(unique(sim$sp_N$cell)), sort(sim$grid_xy$cell))
  expect_equal(sort(unique(sim$sp_N$age)), sim$ages)
  expect_equal(sort(unique(sim$sp_N$year)), sim$years)
  expect_equal(nrow(sim$sp_N), 16 * 3 * 3)
  expect_true(all(is.finite(sim$sp_N$N)))
})

test_that("custom depth closures impose known spatial probabilities", {
  withr::local_seed(438)
  sim <- sim_distribution(small_abundance(), grid = small_grid(),
                           ays_covar = zero_covariance,
                           depth_par = function(data) data$x / 10)
  weights <- exp(sim$grid_xy$x / 10)
  expected <- setNames(weights / sum(weights), sim$grid_xy$cell)
  for (age_value in sim$ages) {
    rows <- sim$sp_N[sim$sp_N$age == age_value & sim$sp_N$year == 2001, ]
    expect_equal(rows$N / sim$N[as.character(age_value), "2001"],
                 unname(expected[as.character(rows$cell)]))
  }
  expect_error(sim_distribution(small_abundance(), grid = small_grid(),
                                ays_covar = zero_covariance,
                                depth_par = function(data) rep(-1e10, nrow(data))),
               "NaN values")
})

test_that("covariance closures honor grouping, dimensions, and reproducibility", {
  withr::local_seed(438)
  x <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  fun <- sim_ays_covar(sd = 0.2, group_ages = 2:3, group_years = 2002:2003)
  result <- fun(x, ages = 1:3, years = 2001:2003, cells = c(2, 4, 6))
  expect_equal(dim(result), c(3L, 3L, 3L))
  expect_equal(dimnames(result), list(age = as.character(1:3),
                                     year = as.character(2001:2003),
                                     cell = c("2", "4", "6")))
  expect_equal(result[2, , ], result[3, , ])
  expect_equal(result[, 2, ], result[, 3, ])
  set.seed(438)
  expect_equal(result, fun(x, ages = 1:3, years = 2001:2003, cells = c(2, 4, 6)))
  expect_error(sim_ays_covar(sd = 1:2)(x, 1:3, 1:3, 1:3), "number of sd")
  expect_error(sim_ays_covar(phi_age = 1:2)(x, 1:3, 1:3, 1:3), "phi_age")
  expect_error(sim_ays_covar(phi_year = 1:2)(x, 1:3, 1:3, 1:3), "phi_year")
})

test_that("optional SPDE covariance has the expected shape", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("INLAspacetime")
  withr::local_seed(438)
  sim <- sim_distribution(small_abundance(), grid = survey_grid,
                           ays_covar = sim_ays_covar_spde(mesh = survey_lite_mesh$mesh,
                                                         model = "spde", sd = 0.1))
  totals <- with(sim$sp_N, tapply(N, list(age, year), sum))
  expect_equal(unname(totals), unname(sim$N), tolerance = 1e-8)
})
