test_that("each fishery stock and ledger use only the corresponding annual population", {
  withr::local_seed(701)
  population <- sim_abundance(ages = 1:3, years = 2001:2003,
    R = sim_R(log_mean = log(c(2000, 1000, 400)), log_sd = 0),
    Z = sim_Z(log_sd = 0)) |>
    sim_distribution(grid = small_grid(), ays_covar = zero_covariance,
      depth_par = function(data) ifelse(data$year == 2002, 0.2 * data$cell, -0.2 * data$cell)) |>
    round_sim()
  saved <- data.table::copy(population)
  key <- .growth_key(population)
  for (year_value in population$years) {
    stock <- .fishery_stock(population, year_value, key)
    actual <- stock[, list(N = sum(N)), by = c("cell", "age")]
    expected <- population$sp_N[population$sp_N$year == year_value, c("cell", "age", "N")]
    data.table::setorderv(actual, c("cell", "age"))
    data.table::setorderv(expected, c("cell", "age"))
    expect_equal(as.data.frame(actual), as.data.frame(expected))
    expect_equal(sum(stock$N), sum(population$N[, as.character(year_value)]))
  }
  specification <- sim_fishery(fishery_share = 0.2, catch_model = sim_capture(1),
                               selectivity = function(age) 1, max_operations = 2)
  result <- sim_fisheries(population, first = specification, second = specification,
                           n_sims = 2)
  ledger <- result$fishery$mortality_ledger
  actual <- ledger[, list(starting = sum(starting), deaths = sum(fishery_deaths + other_deaths),
                           survivors = sum(survivors)), by = c("sim", "year", "cell", "age")]
  for (replicate in 1:2) for (year_value in population$years) {
    rows <- actual[actual$sim == replicate & actual$year == year_value, ]
    expected <- population$sp_N[population$sp_N$year == year_value, c("cell", "age", "N")]
    data.table::setorderv(rows, c("cell", "age"))
    data.table::setorderv(expected, c("cell", "age"))
    expect_equal(rows$starting, expected$N)
    expect_equal(rows$starting, rows$deaths + rows$survivors)
    expect_equal(sum(rows$deaths),
                 round(sum(population$N[, as.character(year_value)] *
                   -expm1(-population$Z[, as.character(year_value)]))))
  }
  expect_equal(population, saved)
})
test_that("fishery location choice follows the current year's shifted stock", {
  withr::local_seed(702)
  population <- fishery_population(100, years = 2001:2002)
  population$sp_N$N[population$sp_N$year == 2002] <- 0
  population$sp_N$N[population$sp_N$year == 2002 & population$sp_N$cell == 2] <- 100
  population <- round_sim(population)
  # A zero-catch callback preserves stocks; no exploration or map error means
  # only the occupied cell can be selected in each year's single operation.
  result <- sim_fishery(population, fishery_share = 0.2,
    catch_model = function(available, operation) numeric(nrow(available)),
    selectivity = function(age) 1, exploration = 0,
    max_operations = 1, sampling_coverage = 0)
  operations <- result$fishery$setdet[order(result$fishery$setdet$year), ]
  expect_equal(operations$year, 2001:2002)
  expect_equal(operations$cell, c(1, 2))
  expect_equal(operations$n, c(0, 0))
})

