biology_population <- function(maturity_by = "length") {
  key <- matrix(c(0.75, 0.25, 0.5, 0.5), 2,
                dimnames = list(length = c(10, 20), age = 1:2))
  sim_abundance(ages = 1:2, years = 2001:2002,
    R = function(years) c(100, 200), N0 = function(R0, Z0, ages) c(100, 40),
    Z = function(ages, years) matrix(log(2), 2, 2, dimnames = list(age = ages, year = years)),
    growth = function(age, length_age_key = FALSE) key,
    weight = function(length) (length / 10)^3,
    maturity = if (maturity_by == "length") function(length) as.numeric(length == 20) else
      function(age) c(0.25, 1), maturity_by = maturity_by)
}

test_that("joint length-age integration produces analytical biomass and mature biomass", {
  sim <- biology_population()
  # W at lengths 10 and 20 is 1 and 8 kg; maturity is 0 and 1.
  expect_equal(sim$B, sim$N * c(2.75, 4.5))
  expect_equal(sim$S, sim$N * c(2, 4))
  expect_equal(sim$B_at_length, sim$N_at_length * c(1, 8))
  expect_equal(sim$S_at_length, sim$N_at_length * c(0, 8))
  expect_equal(colSums(sim$B), colSums(sim$B_at_length))
  expect_equal(colSums(sim$S), colSums(sim$S_at_length))
  expect_equal(dimnames(sim$B), dimnames(sim$N))
  expect_equal(dimnames(sim$S_at_length), dimnames(sim$N_at_length))
  # E[W*M] for age 1 is 2, whereas E[W]*E[M] is only 0.6875.
  expect_false(isTRUE(all.equal(unname(sim$S[1, 1]), 100 * 2.75 * 0.25)))
  expect_true(all(sim$S <= sim$B))
})

test_that("age-based maturity integrates age mixtures into length summaries", {
  sim <- biology_population("age")
  expect_equal(sim$S, sim$B * c(0.25, 1))
  expected <- rbind(0.75 * 0.25 * sim$N[1, ] + 0.5 * sim$N[2, ],
                     8 * (0.25 * 0.25 * sim$N[1, ] + 0.5 * sim$N[2, ]))
  dimnames(expected) <- dimnames(sim$S_at_length)
  expect_equal(sim$S_at_length, expected)
  expect_equal(colSums(sim$S), colSums(sim$S_at_length))
})

test_that("spatial distribution and rounding refresh mass without changing callers", {
  withr::local_seed(402)
  sim <- sim_distribution(biology_population(), grid = small_grid(),
                          ays_covar = zero_covariance, depth_par = function(data) rep(0, nrow(data)))
  saved <- data.table::copy(sim)
  for (metric in c("B", "S")) {
    actual <- tapply(sim$sp_N[[metric]], list(sim$sp_N$age, sim$sp_N$year), sum)
    expect_equal(unname(actual), unname(sim[[metric]]))
  }
  rounded <- round_sim(sim)
  expect_equal(rounded$B, rounded$N * c(2.75, 4.5))
  expect_equal(rounded$S, rounded$N * c(2, 4))
  expect_equal(sum(rounded$sp_N$B), sum(rounded$B))
  expect_equal(sum(rounded$sp_N$S), sum(rounded$S))
  expect_equal(sim, saved)
  expect_equal(dimnames(rounded$B), dimnames(sim$N))
})

test_that("zero abundance and maturity boundaries retain dimensions", {
  sim <- biology_population()
  sim$N[] <- 0
  zero <- .population_biology(sim)
  expect_equal(sum(zero$B), 0)
  expect_equal(sum(zero$S_at_length), 0)
  sim <- biology_population()
  sim$sim_maturity <- function(x) 0
  expect_equal(sum(.population_biology(sim)$S), 0)
  sim$sim_maturity <- function(x) 1
  mature <- .population_biology(sim)
  expect_equal(mature$S, mature$B)
  expect_equal(mature$S_at_length, mature$B_at_length)
  for (ages in list(1, 1:2)) for (years in list(2001, 2001:2002)) {
    one <- sim_abundance(ages = ages, years = years)
    expect_equal(dim(one$B), c(length(ages), length(years)))
    expect_equal(dim(one$S), dim(one$N))
    expect_equal(colSums(one$S), colSums(one$S_at_length))
  }
})

test_that("invalid relationships fail and legacy models are explicit", {
  sim <- biology_population()
  sim$sim_weight <- function(length) -1
  expect_error(.population_biology(sim), "Weight")
  sim$sim_weight <- sim_weight()
  sim$sim_maturity <- function(x) 2
  expect_error(.population_biology(sim), "Maturity")
  sim$sim_maturity <- sim_maturity()
  sim$maturity_by <- "size"
  expect_error(.population_biology(sim), "maturity_by")
  sim$sim_weight <- sim$sim_maturity <- NULL
  expect_identical(.population_biology(sim), sim)
  expect_error(.population_biology(sim, required = TRUE), "migration recipe")
})

test_that("stock_biology summarizes the same moments as B and S", {
  sim <- biology_population()
  out <- stock_biology(sim)
  expect_identical(key(out$at_age), c("year", "age"))
  expect_identical(key(out$at_length), c("year", "length"))
  expect_equal(out$at_age$mean_length, rep(c(12.5, 15), 2))
  expect_equal(out$at_age$mean_weight_kg, rep(c(2.75, 4.5), 2))
  expect_equal(out$at_age$prop_mature, rep(c(0.25, 0.5), 2))
  expect_equal(out$at_age$mean_mature_weight_kg * as.numeric(sim$N), as.numeric(sim$S))
  expect_equal(out$at_length$mean_mature_weight_kg * as.numeric(sim$N_at_length),
               as.numeric(sim$S_at_length))
  sim <- biology_population("age")
  out <- stock_biology(sim)
  expect_equal(out$at_age$prop_mature, rep(c(0.25, 1), 2))
  expect_equal(out$at_length$mean_mature_weight_kg * as.numeric(sim$N_at_length),
               as.numeric(sim$S_at_length))
  sim$N[] <- 0
  expect_equal(stock_biology(sim)$at_age, out$at_age)
  expect_true(all(is.na(stock_biology(sim)$at_length$prop_mature)))
})
