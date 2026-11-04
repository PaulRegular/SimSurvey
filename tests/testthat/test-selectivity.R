test_that("length-based surveys integrate the key and sample selected lengths", {
  withr::local_seed(72)
  sim <- small_population()
  key <- matrix(c(0.75, 0.25, 0.5, 0.5, 0, 1), 2,
                dimnames = list(length = c(10, 20), age = 1:3))
  sim$sim_length <- function(age, length_age_key = FALSE) key
  sim$N_at_length <- convert_N(sim$N, key)
  q <- function(x) as.numeric(x == 20)
  result <- sim_survey(sim, q = q, q_by = "length", trawl_dim = c(5, 5),
                       min_sets = 2, set_den = 0, light = FALSE)
  expect_equal(result$I, result$N * c(0.25, 0.5, 1))
  expect_equal(result$I_at_length[1, ], rep(0, 3), ignore_attr = TRUE)
  expect_equal(result$I_at_length[2, ], colSums(result$I), ignore_attr = TRUE)
  expect_true(all(result$samp$length == 20))
  expect_equal(result$sp_N$I, result$sp_N$N * c(0.25, 0.5, 1)[match(result$sp_N$age, 1:3)])
  expect_error(sim_survey(sim, q_by = "size"), "q_by")
  expect_error(sim_survey(sim, q = function(x) 2, q_by = "length"), "probabilities")
})

test_that("rounding and survey generation leave caller tables unchanged", {
  withr::local_seed(74)
  sim <- small_population()
  saved <- data.table::copy(sim)
  rounded <- round_sim(sim)
  expect_equal(sim, saved)
  result <- sim_survey(sim, trawl_dim = c(1, 5), min_sets = 2, set_den = 0)
  expect_equal(sim, saved)
  expect_null(sim$sp_N$I)
  expect_true(all(rounded$sp_N$N == round(rounded$sp_N$N)))
})

test_that("aggregate sampling conserves counts beyond integer vector sizes", {
  withr::local_seed(32)
  counts <- .multinomial_counts(1e10, c(0.2, 0.3, 0.5))
  expect_equal(sum(counts), 1e10)
  thinned <- .thin_counts(counts, 123)
  expect_equal(sum(thinned), 123)
  expect_true(all(thinned <= counts & thinned >= 0))
  expect_equal(.balanced_round(c(1.4, 1.4, 1.4)), c(2, 1, 1))
})

test_that("the default survey age path is unchanged", {
  withr::local_seed(75)
  sim <- small_population()
  set.seed(76)
  default <- sim_survey(sim, trawl_dim = c(1, 5), min_sets = 2, set_den = 0)
  seed <- .Random.seed
  set.seed(76)
  explicit <- sim_survey(sim, trawl_dim = c(1, 5), min_sets = 2, set_den = 0, q_by = "age")
  expect_equal(default, explicit)
  expect_identical(.Random.seed, seed)
  expect_equal(default$I, default$N * sim_logistic()(sim$ages))
})
