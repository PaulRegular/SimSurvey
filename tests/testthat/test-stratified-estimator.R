test_that("stratified means agree with independently calculated two-stratum totals", {
  data <- data.table::data.table(year = 2001, strat = rep(1:2, each = 2),
                                 strat_area = rep(c(100, 200), each = 2),
                                 tow_area = 10, catch = c(2, 4, 6, 10))
  original <- data.table::copy(data)
  result <- strat_means(data, metric = "catch",
                        strat_groups = c("year", "strat", "strat_area", "tow_area"),
                        survey_groups = "year")
  expected_variance <- ((10 * 8 / 2) * 2 + (20 * 18 / 2) * 8) / 30 ^ 2
  expect_equal(result$mean, 19 / 3)
  expect_equal(result$total, 190)
  expect_equal(result$sd, sqrt(expected_variance))
  expect_equal(result$sets, 4L)
  expect_equal(result$sampling_units, 30)
  expect_true(result$total_lcl < result$total && result$total_ucl > result$total)
  expect_equal(data, original)
})

test_that("zero-variance strata return defined totals and missing confidence limits", {
  data <- data.table::data.table(year = 2001, strat = 1, strat_area = 100,
                                 tow_area = 10, catch = c(0, 0))
  result <- strat_means(data, metric = "catch",
                        strat_groups = c("year", "strat", "strat_area", "tow_area"),
                        survey_groups = "year")
  expect_equal(result$total, 0)
  expect_equal(result$sd, 0)
  expect_true(is.na(result$df))
  expect_true(is.na(result$total_lcl) && is.na(result$total_ucl))
})
