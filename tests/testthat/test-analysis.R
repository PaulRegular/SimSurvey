test_that("error statistics and survey combinations match hand calculations", {
  expect_equal(error_stats(c(-2, 0, 4)),
               c(ME = 2 / 3, MAE = 2, MSE = 20 / 3, RMSE = sqrt(20 / 3)))
  designs <- expand_surveys(set_den = c(0.01, 0.02), lengths_cap = c(5, 10), ages_cap = 1)
  expect_equal(nrow(designs), 4)
  expect_equal(designs$survey, 1:4)
  expect_equal(anyDuplicated(designs[, c("set_den", "lengths_cap", "ages_cap")]), 0L)
})


test_that("full biological sampling reconstructs known lengths and ages", {
  setdet <- data.table::data.table(set = 1:2, sim = 1, year = 2001,
                                   division = 1, strat = 1, n = c(4, 2),
                                   n_measured = c(4, 2), n_aged = c(4, 2))
  samp <- data.table::data.table(set = c(1, 1, 1, 1, 2, 2),
                                 age = c(1, 1, 2, 2, 1, 2),
                                 length = c(10, 10, 20, 20, 10, 20),
                                 measured = TRUE, aged = TRUE)
  result <- strat_data(list(setdet = setdet, samp = samp), length_group = 1)
  expect_equal(result$lf[order(set, length)]$n, c(2, 2, 1, 1))
  expect_equal(result$af[!is.na(age)][order(set, age)]$n, c(2, 2, 1, 1))
  expect_equal(sum(result$lf$n), 6)
  expect_equal(sum(result$af$n), 6)
  setdet$n_aged[2] <- 0
  samp$aged[samp$set == 2] <- FALSE
  expect_warning(strat_data(list(setdet = setdet, samp = samp), alk_scale = "set"),
                 "Age-length keys cannot be constructed")
})

test_that("the complete pipeline retains truth and supports analysis hooks", {
  withr::local_seed(438)
  result <- small_survey(lengths_cap = 1000, ages_cap = 1000,
                         age_sampling = "random") |>
    run_strat() |>
    strat_error()
  expect_equal(nrow(result$total_strat), 3)
  expect_equal(result$I, result$N)
  expect_equal(result$total_strat_error$I, unname(colSums(result$I)))
  expect_true(is.finite(result$total_strat_error_stats["RMSE"]))
  expect_warning(run_strat(result, length_group = 1), "match the length group")

  called <- 0L
  custom_means <- function(...) {
    called <<- called + 1L
    strat_means(...)
  }
  custom <- run_strat(result, strat_means_fun = custom_means)
  expect_equal(called, 3L)
  expect_equal(custom$total_strat, result$total_strat)
})
