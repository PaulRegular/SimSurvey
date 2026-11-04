test_that("survey maps retain all set locations, including zero catches", {
  withr::local_seed(741)
  population <- small_population()
  population$sp_N$N[population$sp_N$cell != 1] <- 0
  population <- round_sim(population)
  sets <- sim_sets(population, cell <= 8, set_den = 0, min_sets = 8,
                    trawl_dim = c(1, 5))
  survey <- sim_survey(population, custom_sets = sets, q = function(age) 1,
                        lengths_cap = 10, ages_cap = 1)
  expect_no_message(plot <- plotly::plotly_build(plot_survey(survey, which_year = 2001)))
  markers <- Filter(function(trace) identical(trace$name, "n") &&
                      identical(trace$mode, "markers"), plot$x$data)[[1]]
  expect_length(markers$x, 8)
  expect_true(all(is.finite(markers$marker$size)))
  expect_equal(markers$marker$sizemin, 5)
  expect_true(any(markers$marker$color == 0))
  expect_equal(.scale_between(rep(0, 3), 2, 600), rep(2, 3))
  expect_equal(.scale_between(rep(5, 3), 2, 600), rep(2, 3))
  expect_equal(.scale_between(numeric()), numeric())
})
