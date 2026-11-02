test_that("survey configurations consume no RNG and require common replicates", {
  withr::local_seed(320)
  rng <- .Random.seed
  config <- sim_survey(q = sim_logistic(x50 = 20), q_by = "length")
  expect_s3_class(config, "sim_survey_spec")
  expect_identical(.Random.seed, rng)
  sim <- small_population()
  expect_error(sim_surveys(sim), "nonempty")
  expect_error(sim_surveys(sim, config), "nonempty")
  expect_error(sim_surveys(sim, a = config, a = config), "unique")
  expect_error(sim_surveys(sim, a = sim_fishery(fishery_share = 0, catch_model = capture_all)), "sim_survey_spec")
  expect_error(sim_surveys(sim, a = sim_survey(n_sims = 2)), "common n_sims")
  sets <- sim_sets(sim, n_sims = 2, min_sets = 2, set_den = 0, trawl_dim = c(1, 5))
  expect_error(sim_surveys(sim, a = sim_survey(custom_sets = sets), n_sims = 1), "replicate identifiers")
  custom <- sim_surveys(sim, a = sim_survey(custom_sets = sets), n_sims = 2)
  expect_equal(sort(custom$survey_programs$a$setdet$set), sort(sets$set))
})

test_that("a single survey program matches the existing standalone interface", {
  withr::local_seed(321)
  sim <- small_population()
  config <- sim_survey(q = function(age) rep(0.8, length(age)), min_sets = 2,
                        trawl_dim = c(1, 5), set_den = 0, lengths_cap = 12,
                        ages_cap = 3, age_sampling = "random", light = FALSE)
  set.seed(322)
  standalone <- sim_survey(sim, q = config$q, min_sets = 2, trawl_dim = c(1, 5),
                           set_den = 0, lengths_cap = 12, ages_cap = 3,
                           age_sampling = "random", light = FALSE)
  rng <- .Random.seed
  set.seed(322)
  batch <- sim_surveys(sim, a = config)
  expect_identical(.Random.seed, rng)
  selected <- select_survey(batch, "a")
  keep <- c("N", "N_at_length", "I", "I_at_length", "sp_N", "setdet", "full_setdet", "samp", "samp_totals")
  expect_equal(selected[keep], standalone[keep])
  expect_null(batch$setdet)
  expect_null(batch$I)
  expect_null(batch$survey_programs$a$N)
  expect_null(batch$survey_programs$a$grid)
  expect_error(select_survey(batch, "missing"), "program")
})

test_that("survey programs are independent of registration and other sampling caps", {
  withr::local_seed(323)
  sim <- small_population()
  a <- sim_survey(q = function(age) 1, trawl_dim = c(1, 5), min_sets = 2,
                   set_den = 0, lengths_cap = 20, ages_cap = 2, light = FALSE)
  b <- a
  b$q <- sim_logistic(x50 = 30)
  b$q_by <- "length"
  set.seed(324)
  full <- sim_surveys(sim, a = a, b = b)
  rng <- .Random.seed
  set.seed(324)
  reordered <- sim_surveys(sim, b = b, a = a)
  expect_equal(full, reordered)
  expect_identical(.Random.seed, rng)
  a$lengths_cap <- 1
  a$ages_cap <- 0
  set.seed(324)
  reduced <- sim_surveys(sim, a = a, b = b)
  for (program in c("a", "b")) {
    expect_equal(full$survey_programs[[program]]$full_setdet,
                 reduced$survey_programs[[program]]$full_setdet)
    expect_equal(full$survey_programs[[program]]$samp$length,
                 reduced$survey_programs[[program]]$samp$length)
  }
  expect_equal(full$survey_programs$b$I_at_length, reduced$survey_programs$b$I_at_length)
})

test_that("survey extraction restores spatial truth and clears cached analyses", {
  withr::local_seed(325)
  sim <- small_population()
  sim$surveys <- expand_surveys(set_den = 0, lengths_cap = 20, ages_cap = 2)
  a <- sim_survey(q = function(age) 1, trawl_dim = c(1, 5), min_sets = 2,
                   set_den = 0, lengths_cap = 100, ages_cap = 4)
  b <- a
  b$q <- function(age) rep(0.5, length(age))
  programs <- sim_surveys(sim, a = a, b = b, n_sims = 2)
  saved <- data.table::copy(programs)
  selected <- select_survey(programs, "a")
  expect_equal(selected$sp_N$I, selected$sp_N$N)
  analysed <- strat_error(run_strat(selected))
  expect_true(nrow(analysed$total_strat) > 0)
  expect_true(nrow(analysed$length_strat_error) > 0)
  expect_s3_class(plot_survey(analysed, which_year = 2001, which_sim = 1), "plotly")
  switched <- select_survey(analysed, "b")
  expect_equal(switched$sp_N$I, switched$sp_N$N * 0.5)
  expect_equal(switched$I, switched$N * 0.5)
  expect_null(switched$total_strat)
  expect_null(switched$age_strat_error)
  expect_null(switched$total_strat_error_stats)
  expect_equal(switched$surveys, sim$surveys)
  expect_equal(programs, saved)
  shuffled <- data.table::copy(programs)
  shuffled$sp_N <- shuffled$sp_N[rev(seq_len(nrow(shuffled$sp_N))), ]
  restored <- select_survey(shuffled, "b")
  expect_equal(restored$sp_N$I, restored$sp_N$N * 0.5)
  switched$setdet$n[] <- 0
  expect_equal(programs, saved)
})

test_that("surveys observe the original population after shared fisheries", {
  withr::local_seed(326)
  sim <- small_population()
  fishery <- sim_fishery(fishery_share = 0.3, catch_model = sim_capture(1),
                          max_operations = 5, lengths_cap = 0)
  fished <- sim_fisheries(sim, first = fishery, second = fishery, n_sims = 2)
  config <- sim_survey(q = function(age) 1, trawl_dim = c(1, 5), min_sets = 2,
                        set_den = 0, lengths_cap = 20, ages_cap = 2)
  set.seed(327)
  surveys <- sim_surveys(fished, a = config)
  set.seed(327)
  original <- sim_surveys(sim, a = config, n_sims = 2)
  expect_equal(surveys$survey_programs, original$survey_programs)
  expect_equal(surveys$fishery, fished$fishery)
  expect_identical(surveys$replicates, 1:2)
  expect_error(sim_surveys(fished, a = config, n_sims = 1), "existing batch")
  expect_error(sim_survey(fished, n_sims = 1), "existing batch")
  expect_equal(sort(unique(sim_survey(fished, trawl_dim = c(1, 5), min_sets = 2,
                  set_den = 0)$setdet$sim)), 1:2)
})

test_that("single-age single-year programs retain matrix dimensions", {
  withr::local_seed(328)
  sim <- sim_abundance(ages = 1, years = 2001,
                        R = sim_R(log_mean = log(1000), log_sd = 0)) |>
    sim_distribution(grid = small_grid(), ays_covar = zero_covariance,
                       depth_par = function(data) rep(0, nrow(data)))
  programs <- sim_surveys(sim,
    a = sim_survey(q = function(age) 1, trawl_dim = c(1, 5), min_sets = 2, set_den = 0),
    b = sim_survey(q = function(length) 0.5, q_by = "length",
                    trawl_dim = c(1, 5), min_sets = 2, set_den = 0))
  a <- select_survey(programs, "a")
  b <- select_survey(programs, "b")
  expect_equal(dim(a$I), c(1L, 1L))
  expect_equal(dimnames(a$I), dimnames(a$N))
  expect_equal(b$I, b$N * 0.5)
  expect_equal(sum(b$I_at_length), sum(b$I))
})
