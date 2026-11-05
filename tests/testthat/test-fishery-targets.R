test_that("target configurations preserve positional calls and the deprecated alias", {
  withr::local_seed(801)
  before <- .Random.seed
  for (type in c("share", "F_full", "catch_n", "catch_kg")) {
    spec <- sim_fishery(target = 0.2, target_type = type, catch_model = capture_all)
    expect_s3_class(spec, "sim_fishery_spec")
    expect_identical(spec$target_type, type)
  }
  expect_identical(.Random.seed, before)
  expect_warning(old <- sim_fishery(fishery_share = 0.2, catch_model = capture_all), "deprecated")
  new <- sim_fishery(target = 0.2, catch_model = capture_all)
  expect_equal(old, new)
  expect_identical(.Random.seed, before)
  expect_error(sim_fishery(target = 0.2, fishery_share = 0.2, catch_model = capture_all), "not both")
  expect_error(sim_fishery(fishery_share = 0.2, target_type = "F_full", catch_model = capture_all), "requires")
  expect_error(sim_fishery(catch_model = capture_all), "Supply target")
  for (type in list("Fbar", c("share", "F_full"), NA_character_, 1)) {
    expect_error(sim_fishery(target = 0.2, target_type = type, catch_model = capture_all), "target_type")
  }
  sim <- fishery_population()
  conflicting <- new
  conflicting$fishery_share <- 0.3
  expect_error(sim_fisheries(sim, net = conflicting), "net.*not both")
  legacy <- new
  legacy$fishery_share <- legacy$target
  legacy$target <- legacy$target_type <- NULL
  set.seed(802)
  a <- sim_fisheries(sim, fishery = legacy)$fishery
  set.seed(802)
  b <- sim_fisheries(sim, fishery = new)$fishery
  expect_equal(a$setdet, b$setdet)
  expect_equal(a$reference, b$reference)
})

test_that("fully selected F preserves selectivity magnitude and analytical references", {
  withr::local_seed(803)
  sim <- fishery_population()
  sim$sim_weight <- function(length) length / 10
  sim <- round_sim(sim)
  result <- sim_fishery(sim, 0.4, capture_all, target_type = "F_full",
                         selectivity = function(age) age / 2, subset_cells = cell == 1)$fishery
  expected_n <- 100 * c(0.2, 0.4) / log(2) * 0.5
  expect_equal(result$reference$F_reference, c(0.2, 0.4))
  expect_equal(result$reference$catch_reference, expected_n)
  expect_equal(result$reference$kg_reference, expected_n * c(1, 2))
  expect_equal(result$annual$F_full_reference, 0.4)
  expect_equal(result$annual$requested_target, 0.4)
  expect_equal(result$annual$target, round(sum(expected_n)))
  expect_equal(result$annual$target_unit, "fish")
  expect_equal(result$annual$n, result$annual$target)
  expect_equal(result$annual$achieved, result$annual$n)
  # F_full need not be below juvenile Z when juvenile selectivity is small.
  sim$Z[, 1] <- c(0.1, 0.4)
  result <- sim_fishery(sim, 0.2, capture_all, target_type = "F_full",
                         selectivity = function(age) c(0.25, 1))$fishery
  expect_equal(result$reference$F_reference, c(0.05, 0.2))
})

test_that("number and kilogram targets imply consistent reference F", {
  withr::local_seed(804)
  sim <- fishery_population()
  sim$sim_weight <- function(length) length / 10
  sim <- round_sim(sim)
  run <- function(target, type) sim_fishery(sim, target, capture_all, target_type = type,
                                            selectivity = function(age) 1,
                                            subset_cells = cell == 1, light = FALSE)$fishery
  numbers <- run(30, "catch_n")
  biomass <- run(45, "catch_kg")
  expect_equal(numbers$reference$F_reference, rep(0.3 * log(2), 2))
  reference_columns <- c("N", "Z", "selectivity", "F_reference", "catch_reference", "kg_reference")
  expect_equal(biomass$reference[, ..reference_columns], numbers$reference[, ..reference_columns])
  expect_equal(numbers$annual$target, 30)
  expect_equal(biomass$annual$target, 45)
  expect_equal(biomass$annual$target_unit, "kg")
  expect_equal(biomass$annual$reference_target_n, 30)
  expect_equal(biomass$annual$reference_target_kg, 45)
  expect_equal(biomass$annual$kg, 45)
  expect_equal(biomass$annual$achieved, biomass$annual$kg)
  expect_equal(biomass$annual$shortfall, 0)
  expect_equal(biomass$annual$status, "target_met")
  expect_true(biomass$setdet$quota_limited)
  expect_equal(sum(biomass$joint_catch$kg), 45)
  expect_true(all(biomass$joint_catch$n == floor(biomass$joint_catch$n)))
})

test_that("biomass targets integrate weight jointly with length selectivity", {
  withr::local_seed(805)
  sim <- fishery_population(100)
  key <- matrix(c(0.75, 0.25), 2, dimnames = list(length = c(10, 20), age = 1))
  sim$sim_length <- function(age, length_age_key = FALSE) key
  sim$sim_weight <- function(length) ifelse(length == 10, 1, 8)
  result <- sim_fishery(sim, 10, capture_all, target_type = "catch_kg",
                         selectivity = function(length) ifelse(length == 20, 0.5, 0),
                         selectivity_by = "length", subset_cells = cell == 1, light = FALSE)$fishery
  # E[w*s] = .25*8*.5 = 1, rather than E[w]*E[s] = 2.75*.125.
  expect_equal(result$annual$F_full_reference, 0.2 * log(2))
  expect_equal(result$reference$selectivity, 0.125)
  expect_equal(result$reference$catch_reference, 1.25)
  expect_equal(result$reference$kg_reference, 10)
  expect_equal(result$annual$kg, 8)
  expect_equal(result$annual$shortfall, 2)
  expect_equal(result$annual$status, "exhausted")
  expect_equal(result$joint_catch$length, 20)
})

test_that("invalid targets report population and mortality limits before RNG draws", {
  withr::local_seed(806)
  sim <- fishery_population(100)
  sim$sim_weight <- function(length) 2
  sim <- round_sim(sim)
  spec <- function(value, type) sim_fishery(target = value, target_type = type,
                                             catch_model = capture_all, selectivity = function(age) 1)
  before <- .Random.seed
  expect_error(sim_fisheries(sim, trawl = spec(101, "catch_n")), "trawl.*2001.*101.*N 100")
  expect_error(sim_fisheries(sim, trawl = spec(201, "catch_kg")), "trawl.*2001.*201.*B 200")
  expect_error(sim_fisheries(sim, trawl = spec(60, "catch_n")), "trawl.*2001.*60.*mortality-compatible limit 50")
  expect_error(sim_fisheries(sim, trawl = spec(120, "catch_kg")), "trawl.*2001.*120.*mortality-compatible limit 100")
  expect_error(sim_fisheries(sim, trawl = spec(0.8, "F_full")), "trawl.*2001.*age 1.*exceeds Z.*maximum F_full")
  expect_error(sim_fisheries(sim, trawl = spec(0.4, "F_full"), net = spec(30, "catch_n")),
               "Combined fisheries.*net.*trawl.*2001.*age 1.*exceeds Z")
  for (value in list(-1, Inf, NA_real_, numeric(), "10")) {
    expect_error(sim_fisheries(sim, net = spec(value, "catch_n")), "net.*target")
  }
  expect_identical(.Random.seed, before)
  sim <- fishery_population(100, years = 2001:2002)
  before <- .Random.seed
  expect_error(sim_fisheries(sim, net = spec(c(0.1, 0.8), "F_full")), "net.*2002.*age 1.*exceeds Z")
  expect_identical(.Random.seed, before)
})

test_that("targets align years and cover zero populations, selection and mortality", {
  withr::local_seed(807)
  sim <- fishery_population(100, years = 2001:2002)
  result <- sim_fishery(sim, c(`2002` = 0.4, `2001` = 0.1), capture_all, target_type = "F_full",
                         selectivity = function(age) 1)$fishery
  expect_equal(result$reference$F_reference, c(0.1, 0.4))
  expect_equal(result$annual$requested_target, c(0.1, 0.4))
  expect_error(sim_fishery(sim, c(`2002` = 1, `2003` = 2), capture_all, target_type = "catch_n"), "names must match")
  expect_error(sim_fishery(sim, 1:3, capture_all, target_type = "catch_kg"), "target")
  for (type in c("share", "F_full", "catch_n", "catch_kg")) {
    empty <- sim_fishery(fishery_population(0), 0, capture_all, target_type = type)$fishery
    expect_equal(empty$annual$achieved, 0)
    expect_equal(empty$annual$status, "target_met")
    expect_true(is.na(empty$reference$F_realized))
  }
  sim <- fishery_population(100)
  sim$Z[] <- 0
  for (type in c("F_full", "catch_n", "catch_kg")) {
    expect_error(sim_fishery(sim, 0.1, capture_all, target_type = type), "exceeds")
  }
  expect_equal(sim_fishery(sim, 1, capture_all)$fishery$annual$n, 0)
  sim <- fishery_population(100)
  expect_error(sim_fishery(sim, 1, capture_all, target_type = "catch_n",
                            selectivity = function(age) 0), "limit 0")
  expect_equal(sim_fishery(sim, 1, capture_all, target_type = "F_full",
                            selectivity = function(age) 0)$fishery$annual$n, 0)
  sim$sim_weight <- function(length) 0
  expect_error(sim_fishery(sim, 1, capture_all, target_type = "catch_kg"), "B 0")
  sim <- fishery_population(100)
  sim$Z[] <- 1e-14
  tiny <- sim_fishery(sim, 5e-15, capture_all, target_type = "F_full",
                       selectivity = function(age) 1)$fishery
  expect_equal(tiny$reference$catch_reference, 100 * 0.5 * -expm1(-1e-14))
  expect_error(sim_fishery(sim, 2e-14, capture_all, target_type = "F_full",
                            selectivity = function(age) 1), "exceeds Z")
  expect_error(sim_fishery(sim, 2e-12, capture_all, target_type = "catch_n",
                            selectivity = function(age) 1), "mortality-compatible limit")
  sim$sim_weight <- function(length) 1e-20
  expect_error(sim_fishery(sim, 2e-18, capture_all, target_type = "catch_kg"), "exceeds starting B")
})

test_that("all target modes share stocks and are invariant to registration order", {
  withr::local_seed(808)
  sim <- fishery_population()
  sim$sim_weight <- function(length) 1
  configs <- lapply(c("share", "F_full", "catch_n", "catch_kg"), function(type) {
    sim_fishery(target = if (type %in% c("share", "F_full")) 0.1 else 10,
                  target_type = type, catch_model = capture_all, selectivity = function(age) 1,
                  subset_cells = cell == 1, light = FALSE)
  })
  names(configs) <- c("share", "full", "numbers", "biomass")
  set.seed(809)
  first <- do.call(sim_fisheries, c(list(sim = sim, n_sims = 2), configs))$fishery
  rng <- .Random.seed
  set.seed(809)
  second <- do.call(sim_fisheries, c(list(sim = sim, n_sims = 2), rev(configs)))$fishery
  expect_equal(first, second)
  expect_identical(.Random.seed, rng)
  ledger <- first$mortality_ledger
  expect_equal(ledger$starting, ledger$fishery_deaths + ledger$other_deaths + ledger$survivors)
  expect_true(all(ledger$survivors >= 0))
  expect_true(all(first$annual$status == "target_met"))
  expect_true(all(first$annual$shortfall == 0))
  for (type in c("F_full", "catch_n", "catch_kg")) {
    target <- if (type == "F_full") 0.1 else 10
    sim$replicates <- 1:2
    set.seed(810)
    direct <- sim_fishery(sim, target, capture_all, target_type = type)
    rng <- .Random.seed
    set.seed(810)
    batch <- sim_fisheries(sim, fishery = sim_fishery(target = target, catch_model = capture_all,
                                                     target_type = type))
    expect_equal(direct, batch)
    expect_identical(.Random.seed, rng)
  }
})

test_that("kilogram quotas stay aggregate, bounded and independent of monitoring", {
  withr::local_seed(811)
  counts <- c(1e10, 2e10, 3e10)
  weights <- c(0, 1.3, 4.7)
  for (limit in c(0, 0.2, 10.5, 1e9, sum(counts * weights))) {
    caught <- .thin_biomass(counts, weights, limit)
    expect_true(all(caught >= 0 & caught <= counts & caught == floor(caught)))
    expect_lte(sum(caught * weights), limit)
    expect_equal(caught[1], counts[1])
  }
  sim <- fishery_population(1e10)
  sim$sim_weight <- function(length) 3
  result <- sim_fishery(sim, 1e9, capture_all, target_type = "catch_kg",
                         selectivity = function(age) 1, subset_cells = cell == 1,
                         sampling_coverage = 0, light = FALSE)$fishery
  expect_equal(result$annual$n, floor(1e9 / 3))
  expect_equal(result$annual$shortfall, 1)
  expect_equal(nrow(result$joint_catch), 1)
  expect_equal(nrow(result$samp), 0)
  sim <- fishery_population()
  saved <- data.table::copy(sim)
  run <- function(...) sim_fishery(sim, 0.01, sim_capture(1), target_type = "catch_kg",
                                    max_operations = 5, light = FALSE, ...)$fishery
  set.seed(812)
  a <- run()
  set.seed(812)
  b <- run(sampling_coverage = 0, reporting_cv = 1)
  expect_equal(a$joint_catch, b$joint_catch)
  expect_equal(a$reference, b$reference)
  expect_equal(a$mortality_ledger, b$mortality_ledger)
  expect_equal(sim, saved)
})

test_that("biomass target plots use kilogram targets and shortfall units", {
  withr::local_seed(813)
  sim <- fishery_population(100)
  sim$sim_weight <- function(length) 2
  result <- sim_fishery(sim, 21, capture_all, target_type = "catch_kg",
                         selectivity = function(age) 1, subset_cells = cell == 1)
  kg <- plotly::plotly_build(plot_fishery_trend(result, "fishery", "kg"))
  expect_equal(as.numeric(kg$x$data[[3]]$y), 21)
  n <- plotly::plotly_build(plot_fishery_trend(result, "fishery", "n"))
  expect_length(n$x$data, 2)
  shortfall <- plotly::plotly_build(plot_fishery_trend(result, "fishery", "shortfall"))
  expect_match(paste(unlist(shortfall$x$layout$yaxis$title), collapse = " "), "kg")
  expect_match(shortfall$x$data[[1]]$text, "Shortfall \\(kg\\)")
})

test_that("kilogram targets retain unsuccessful operations and repeated-cell depletion", {
  withr::local_seed(814)
  sim <- fishery_population(100)
  sim$sim_weight <- function(length) 2
  capture <- function(available, operation) {
    if (operation$fishery_operation == 1) return(available$N * 0)
    .thin_counts(available$N, 4)
  }
  result <- sim_fishery(sim, 21, capture, target_type = "catch_kg",
                         selectivity = function(age) 1, subset_cells = cell == 1,
                         light = FALSE)$fishery
  expect_equal(result$setdet$n, c(0, 4, 4, 2))
  expect_equal(result$setdet$kg, c(0, 8, 8, 4))
  expect_equal(result$setdet$quota_limited, c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(result$annual$effort, 4)
  expect_equal(result$annual$shortfall, 1)
  expect_equal(result$annual$status, "exhausted")
  ledger <- result$mortality_ledger
  expect_equal(ledger$starting, ledger$fishery_deaths + ledger$other_deaths + ledger$survivors)
  expect_equal(sum(ledger$fishery_deaths), 10)
  expect_equal(sum(ledger$other_deaths), 40)
  expect_equal(sum(ledger$survivors), 50)
  limited <- sim_fishery(sim, 21, capture, target_type = "catch_kg",
                          selectivity = function(age) 1, subset_cells = cell == 1,
                          max_operations = 1)$fishery
  expect_equal(limited$annual$status, "operation_limit")
  expect_equal(limited$annual$shortfall, 21)
  expect_equal(limited$annual$kg, 0)
})
