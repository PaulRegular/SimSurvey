test_that("parallel survey batches have unique set and replicate identifiers", {
  withr::local_seed(438)
  local_worker_library()
  for (cores in c(1, 2)) {
    result <- sim_survey_parallel(small_population(), n_sims = 2, n_loops = 2,
                                   cores = cores, quiet = TRUE,
                                   q = function(x) rep(1, length(x)),
                                   trawl_dim = c(1, 5), min_sets = 2, set_den = 0,
                                   lengths_cap = 10, ages_cap = 1)
    expect_equal(sort(unique(result$setdet$sim)), 1:4)
    expect_equal(anyDuplicated(result$setdet$set), 0L)
    expect_true(all(result$samp$set %in% result$setdet$set))
    expect_equal(nrow(result$samp_totals), 4 * 3)
    expect_equal(sum(result$samp_totals$n_caught), sum(result$setdet$n))
  }
})

test_that("survey scenarios export and resume in temporary directories", {
  withr::local_seed(438)
  local_worker_library()
  export_dir <- withr::local_tempdir()
  designs <- expand_surveys(set_den = 0, lengths_cap = c(30, 40), ages_cap = 2)
  result <- test_surveys(small_population(), surveys = designs,
                         n_sims = 1, n_loops = 1, cores = 1,
                         export_dir = export_dir, progress = FALSE,
                         q = function(x) rep(1, length(x)), trawl_dim = c(1, 5),
                         age_sampling = "random")
  expect_true(all(file.exists(file.path(export_dir,
                                        c("test_inputs.RData", "complete.RData", "test_output.RData")))))
  expect_equal(sort(unique(result$samp_totals$survey)), 1:2)
  expected_length_truth <- convert_N(result$I,
                                      result$sim_length(result$ages, length_age_key = TRUE))
  expect_equal(result$I_at_length, expected_length_truth)
  resumed <- resume_test(export_dir)
  expect_equal(resumed$I_at_length, expected_length_truth)
  expect_equal(resumed$samp_totals, result$samp_totals)
  expect_equal(resumed$total_strat_error, result$total_strat_error)

  state <- new.env(parent = emptyenv())
  load(file.path(export_dir, "complete.RData"), envir = state)
  expect_true(all(state$complete))
  state$complete[2] <- FALSE
  save(list = "complete", file = file.path(export_dir, "complete.RData"), envir = state)
  partial <- resume_test(export_dir)
  expect_equal(partial$samp_totals[partial$samp_totals$survey == 1, ],
               result$samp_totals[result$samp_totals$survey == 1, ])
  expect_equal(sort(unique(partial$samp_totals$survey)), 1:2)
  expect_equal(nrow(partial$samp_totals), nrow(result$samp_totals))
})

test_that("survey programs coexist with scenario designs and length-selective wrappers", {
  withr::local_seed(330)
  local_worker_library()
  sim <- small_population()
  q <- function(length) as.numeric(length >= 30)
  config <- sim_survey(q = q, q_by = "length", trawl_dim = c(1, 5),
                        min_sets = 2, set_den = 0, lengths_cap = 30, ages_cap = 2)
  programs <- sim_surveys(sim, a = config)
  designs <- expand_surveys(set_den = 0, lengths_cap = 30, ages_cap = 2)
  result <- test_surveys(programs, surveys = designs, n_sims = 1, n_loops = 1,
                         cores = 1, progress = FALSE, q = q, q_by = "length",
                         trawl_dim = c(1, 5))
  expect_equal(result$surveys, designs)
  expect_equal(result$survey_programs, programs$survey_programs)
  view <- select_survey(result, "a")
  expect_null(view$total_strat_error)
  expect_equal(view$surveys, designs)
  expect_equal(view$I_at_length, programs$survey_programs$a$I_at_length)
  parallel <- sim_survey_parallel(sim, n_sims = 1, n_loops = 1, cores = 1, quiet = TRUE,
                                   q = q, q_by = "length", trawl_dim = c(1, 5),
                                   min_sets = 2, set_den = 0, lengths_cap = 30, ages_cap = 2)
  expect_equal(parallel$I, programs$survey_programs$a$I)
  expect_equal(parallel$I_at_length, programs$survey_programs$a$I_at_length)
})

test_that("maturity and refreshed biomass propagate through parallel and scenario wrappers", {
  withr::local_seed(430)
  local_worker_library()
  sim <- small_population()
  expected <- round_sim(sim)
  parallel <- sim_survey_parallel(sim, n_sims = 1, n_loops = 1, cores = 1, quiet = TRUE,
    q = function(age) 1, trawl_dim = c(1, 5), min_sets = 2, set_den = 0,
    lengths_cap = 20, ages_cap = 2, maturity_sampling = "measured")
  expect_equal(parallel$B, expected$B)
  expect_equal(parallel$S, expected$S)
  expect_equal(sum(parallel$sp_N$B), sum(parallel$B))
  expect_equal(parallel$setdet$n_maturity, parallel$setdet$n_measured)
  expect_equal(parallel$samp$maturity_observed, parallel$samp$measured)
  designs <- expand_surveys(set_den = 0, lengths_cap = 30, ages_cap = 2)
  export_dir <- withr::local_tempdir()
  scenarios <- test_surveys(sim, surveys = designs, n_sims = 1, n_loops = 1,
    cores = 1, progress = FALSE, export_dir = export_dir, q = function(age) 1,
    trawl_dim = c(1, 5), maturity_sampling = "measured")
  expect_equal(scenarios$B, expected$B)
  expect_equal(scenarios$S, expected$S)
  expect_equal(scenarios$samp_totals$n_maturity, scenarios$samp_totals$n_measured)
  expect_equal(scenarios$setdet$n_maturity, scenarios$setdet$n_measured)
  resumed <- resume_test(export_dir)
  expect_equal(resumed$S, scenarios$S)
  columns <- c("set", "id", "age", "length", "mature")
  expect_equal(resumed$sample_truth[, ..columns], resumed$samp[resumed$samp$measured, ..columns])
  expect_equal(resumed$sample_truth$true_mature, resumed$sample_truth$mature)
  expect_equal(resumed$samp_totals, scenarios$samp_totals)
})
