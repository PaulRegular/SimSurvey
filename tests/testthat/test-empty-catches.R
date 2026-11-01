test_that("catchability extremes retain zero catches and bound full catches", {
  withr::local_seed(438)
  expect_equal(group_lengths(numeric(), 3), numeric())
  sim <- small_population()
  sets <- sim_sets(sim, set_den = 0, min_sets = 2, trawl_dim = c(5, 5))
  zero <- sim_survey(sim, custom_sets = sets,
                     q = function(x) rep(0, length(x)), lengths_cap = 5, ages_cap = 1)
  expect_equal(nrow(zero$setdet), nrow(sets))
  expect_true(all(zero$setdet$n == 0))
  expect_true(all(zero$setdet$n_measured == 0 & zero$setdet$n_aged == 0))
  expect_equal(nrow(zero$samp), 0)
  full <- sim_survey(sim, custom_sets = sets,
                     q = function(x) rep(1, length(x)), binom_error = FALSE,
                     lengths_cap = 5, ages_cap = 1)
  expect_equal(full$setdet$n, full$setdet$N)
  expect_equal(full$I, full$N)
})

test_that("zero biological sampling retains positive catches without samples", {
  withr::local_seed(438)
  result <- small_survey(lengths_cap = 0, ages_cap = 0, age_sampling = "random")
  expect_true(any(result$setdet$n > 0))
  expect_true(all(result$setdet$n_measured == 0 & result$setdet$n_aged == 0))
  expect_false(any(result$samp$measured | result$samp$aged))
})
