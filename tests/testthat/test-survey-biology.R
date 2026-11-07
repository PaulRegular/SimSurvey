survey_biology_fixture <- function(protocol = "measured") {
  setdet <- data.table::data.table(
    set = 1:4, sim = 1, year = 2001, division = 1,
    strat = c(1, 1, 2, 2), strat_area = c(40, 40, 45, 45),
    tow_area = c(2, 2, 5, 5), n = c(12, 6, 8, 0))
  samp <- data.table::data.table(
    set = c(1, 1, 1, 2, 2, 3, 3, 1),
    length = c(10, 10, 20, 10, 20, 10, 20, 999),
    age = c(1, 99, 2, 2, 99, 1, 2, 99),
    weight = c(2, 4, 8, 6, 10, 14, 16, Inf),
    mature = c(0, 1, 1, 1, 0, 0, 0, 99),
    measured = c(rep(TRUE, 7), FALSE),
    aged = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE),
    age_probability = c(rep(2 / 3, 3), 0.5, 0.5, 1, 1, NA_real_))
  samp$maturity_observed <- switch(protocol, measured = samp$measured, aged = samp$aged,
                                   none = rep(FALSE, nrow(samp)))
  samp$maturity_probability <- switch(protocol, measured = ifelse(samp$measured, 1, NA_real_),
    aged = samp$age_probability, none = ifelse(samp$measured, 0, NA_real_))
  list(setdet = setdet, samp = samp, ages = 1:3, lengths = c(10, 20, 30),
       years = 2001, sim_length = sim_vonB(length_group = 1))
}

test_that("survey trait ratios use swept-area and conditional sampling weights", {
  measured <- run_strat(survey_biology_fixture("measured"))
  aged <- run_strat(survey_biology_fixture("aged"))
  none <- run_strat(survey_biology_fixture("none"))

  ## The tow expansion factors are 10, 10 and 4.5. A measured fish represents
  ## 4, 3 and 4 caught fish in those sets. Aging further expands by 1.5, 2 and 1.
  expect_equal(measured$length_strat$mean_weight_kg, c(672 / 128, 908 / 88, NA))
  expect_equal(measured$length_strat$prop_mature, c(70 / 128, 40 / 88, NA))
  expect_equal(aged$length_strat$prop_mature, c(60 / 138, 60 / 78, NA))
  expect_equal(measured$age_strat[!is.na(age)]$mean_weight_kg,
               c((60 * 2 + 18 * 14) / 78, (60 * 8 + 60 * 6 + 18 * 16) / 138, NA))
  expect_equal(measured$age_strat[!is.na(age)]$prop_mature, c(0, 120 / 138, NA))
  expect_equal(aged$age_strat$prop_mature, measured$age_strat$prop_mature)
  expect_equal(aged$length_strat$mean_weight_kg, measured$length_strat$mean_weight_kg)
  expect_equal(none$length_strat$mean_weight_kg, measured$length_strat$mean_weight_kg)
  expect_true(all(is.na(none$length_strat$prop_mature)))
  expect_identical(none$length_strat$maturity_status,
                   c("no_maturity_support", "no_maturity_support", "no_observations"))
  expect_equal(measured$age_strat[age == 1]$maturity_status, "ok")
  expect_equal(measured$age_strat[is.na(age)]$maturity_status, "unassigned_age")
  expect_equal(measured$length_strat[length == 30]$weight_status, "no_observations")
  expect_false("mean_weight_kg" %in% names(measured$total_strat))
  expect_equal(aged$age_strat$total, measured$age_strat$total)
  expect_equal(none$length_strat$total, measured$length_strat$total)

  altered <- survey_biology_fixture("aged")
  altered$samp$age[!altered$samp$aged] <- -999
  altered$samp$mature[!altered$samp$maturity_observed] <- -999
  altered$samp$weight[!altered$samp$measured] <- -999
  altered$setdet$kg <- -999
  altered <- run_strat(altered)
  expect_equal(altered$age_strat, aged$age_strat)
  expect_equal(altered$length_strat, aged$length_strat)
})

test_that("missing biological design or observations do not change count estimates", {
  survey <- survey_biology_fixture()
  result <- run_strat(survey)
  legacy <- data.table::copy(survey)
  legacy$samp$age_probability <- legacy$samp$maturity_probability <- NULL
  estimated <- run_strat(legacy)
  expect_equal(estimated$length_strat$total, result$length_strat$total)
  expect_equal(estimated$age_strat$total, result$age_strat$total)
  expect_equal(estimated$length_strat$mean_weight_kg, result$length_strat$mean_weight_kg)
  expect_true(all(is.na(estimated$age_strat$mean_weight_kg)))
  expect_true(all(estimated$age_strat[!is.na(age)]$weight_status == "age_design_unknown"))
  expect_true(all(is.na(estimated$length_strat$prop_mature)))

  missing <- data.table::copy(survey)
  missing$samp$weight[1] <- NA_real_
  missing$samp$maturity_observed[2] <- FALSE
  missing <- run_strat(missing)
  expect_identical(missing$length_strat$weight_status, c("missing_weight", "ok", "no_observations"))
  expect_identical(missing$length_strat$maturity_status,
                   c("missing_maturity", "ok", "no_observations"))
  expect_equal(missing$age_strat[age == 1]$maturity_status, "ok")

  aged_only <- survey_biology_fixture("aged")
  aged_only$samp$maturity_observed[1] <- FALSE
  aged_only <- run_strat(aged_only)
  expect_equal(aged_only$length_strat[length == 10]$maturity_status, "missing_maturity")
  expect_equal(aged_only$age_strat[age == 1]$maturity_status, "missing_maturity")

  custom <- run_strat(survey, strat_means_fun = function(...) strat_means(...))
  expect_equal(custom$age_strat$total, result$age_strat$total)
  expect_true(all(custom$age_strat[!is.na(age)]$weight_status == "custom_estimator"))
  expect_true(all(is.na(custom$length_strat$prop_mature)))
})

test_that("survey sampling records both conditional probabilities", {
  withr::local_seed(821)
  for (design in c("random", "stratified")) {
    survey <- small_survey(age_sampling = design, age_space_group = "division",
                           age_length_group = 4, lengths_cap = 16, ages_cap = 2,
                           maturity_sampling = "aged", light = FALSE)
    expect_true(all(is.na(survey$samp$age_probability[!survey$samp$measured])))
    expect_true(all(is.na(survey$samp$maturity_probability[!survey$samp$measured])))
    measured <- merge(survey$samp[survey$samp$measured, ],
                      survey$setdet[, c("set", "sim", "year", "division")], by = "set")
    measured$bin <- group_lengths(measured$length, 4)
    groups <- if (design == "random") "set" else c("sim", "year", "division", "bin")
    expected <- measured[, list(expected = rep(min(.N, 2) / .N, .N), actual = age_probability),
                         by = groups]
    expect_equal(expected$actual, expected$expected)
    expect_equal(measured$maturity_probability, measured$age_probability)
    expect_true(all(measured$age_probability > 0))
    expect_true(any(measured$age_probability < 1))
  }
  zero <- small_survey(age_sampling = "random", ages_cap = 0, maturity_sampling = "measured")
  expect_true(all(zero$samp$age_probability[zero$samp$measured] == 0))
  expect_true(all(zero$samp$maturity_probability[zero$samp$measured] == 1))
  zero <- run_strat(zero)
  expect_true(all(zero$age_strat[!is.na(age)]$weight_status == "no_aging_support"))
  expect_true(any(zero$length_strat$maturity_status == "ok"))
})
