test_that("rounding produces internally consistent population summaries", {
  withr::local_seed(438)
  sim <- round_sim(small_population())
  expect_equal(sim$sp_N$N, round(sim$sp_N$N))
  totals <- with(sim$sp_N, tapply(N, list(age, year), sum))
  expect_equal(unname(totals), unname(sim$N))
  expect_equal(colSums(sim$N_at_length), colSums(sim$N))
  expect_equal(unname(sim$N0), unname(sim$N[, 1]))
  expect_equal(unname(sim$R), unname(sim$N[1, ]))
})

test_that("set allocation honors density, minimums, subsets, and repeats", {
  withr::local_seed(438)
  sim <- small_population()
  sets <- sim_sets(sim, n_sims = 2, set_den = 0, min_sets = 2)
  expect_equal(nrow(sets), 2 * 3 * 2)
  expect_equal(anyDuplicated(sets$set), 0L)
  expect_true(all(sets$cell_sets == 1))
  expect_true(all(sets$strat_area == 400))
  expect_equal(sim_sets(sim, set_den = 0.01, min_sets = 0)$strat_sets, rep(4, 12))
  subset <- sim_sets(sim, cell <= 4, set_den = 0, min_sets = 2)
  expect_true(all(subset$cell <= 4))
  repeated <- sim_sets(sim, cell == 1, set_den = 0, min_sets = 3,
                       resample_cells = TRUE)
  expect_true(all(repeated$cell_sets == 3))
  expect_error(sim_sets(sim, cell == 1, min_sets = 3, set_den = 0), "sample")
})


test_that("custom sets and random biological samples honor caps", {
  withr::local_seed(438)
  sim <- small_population()
  sets <- sim_sets(sim, n_sims = 2, min_sets = 2, set_den = 0,
                   trawl_dim = c(1, 5))
  result <- sim_survey(sim, n_sims = 2, custom_sets = sets,
                       q = function(x) rep(1, length(x)), lengths_cap = 12,
                       ages_cap = 3, age_sampling = "random", light = FALSE)
  expect_equal(sort(result$setdet$set), sort(sets$set))
  expect_true(all(result$setdet$n_measured <= 12))
  expect_true(all(result$setdet$n_aged <= 3))
  expect_true(all(!result$samp$aged | result$samp$measured))
  expect_equal(nrow(result$samp), sum(result$setdet$n))
  expect_true(all(result$full_setdet$n <= result$full_setdet$N))
  expect_equal(sum(result$samp$measured), sum(result$setdet$n_measured))
  expect_equal(sum(result$samp$aged), sum(result$setdet$n_aged))
  duplicated <- data.table::copy(sets)
  duplicated$set[2] <- duplicated$set[1]
  expect_error(sim_survey(sim, custom_sets = duplicated), "unique number")
})

test_that("stratified age samples respect each configured sampling scale", {
  withr::local_seed(438)
  for (scale in c("division", "strat", "set")) {
    result <- small_survey(lengths_cap = 30, ages_cap = 1,
                            age_sampling = "stratified", age_space_group = scale)
    aged <- merge(result$samp[result$samp$aged, ], result$setdet,
                   by = "set", suffixes = c("", "_set"))
    aged$length_group <- group_lengths(aged$length, 1)
    groups <- aged[, .N, by = c("sim", "year", scale, "length_group")]
    expect_true(all(groups$N <= 1))
    expect_true(all(!result$samp$aged | result$samp$measured))
  }
  expect_error(small_survey(age_sampling = "other"), "age_sampling")
  expect_error(small_survey(age_space_group = "other"), "age_space_group")
  expect_error(small_survey(age_sampling = "random", lengths_cap = 2, ages_cap = 3),
               "cannot exceed")
})

test_that("repeated sets share a bounded cell catch", {
  withr::local_seed(438)
  sim <- small_population()
  sets <- sim_sets(sim, cell == 1, min_sets = 3, set_den = 0,
                   resample_cells = TRUE, trawl_dim = c(5, 5))
  result <- sim_survey(sim, custom_sets = sets, q = function(x) rep(1, length(x)),
                       binom_error = FALSE, lengths_cap = 10, ages_cap = 1,
                       light = FALSE)
  totals <- result$full_setdet[, list(caught = sum(n), available = unique(N)),
                               by = c("year", "cell", "age")]
  expect_true(all(totals$caught <= totals$available + 1))
})
