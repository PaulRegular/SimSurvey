test_that("batch estimates and both truth scopes match single-program calls", {
  withr::local_seed(1204)
  population <- small_population()
  population$sp_N$N <- population$sp_N$N * (1 + population$sp_N$cell %% 3)
  population <- round_sim(population)
  west <- sim_sets(population, (sim == 1 & x <= 5) | (sim == 2 & x <= 0 & year != 2002),
                   n_sims = 2, min_sets = 1, set_den = 0, trawl_dim = c(1, 5))
  east <- sim_sets(population, x >= -5, n_sims = 2, min_sets = 1,
                   set_den = 0, trawl_dim = c(1, 5))
  age_q <- function(age) age / 4
  length_q <- function(length) ifelse(length <= 30, 0.25, 0.8)
  batch <- population |>
    sim_surveys(west = sim_survey(q = age_q, custom_sets = west),
                  east = sim_survey(q = length_q, q_by = "length", custom_sets = east),
                  n_sims = 2)
  original <- data.table::copy(batch)
  estimates <- run_strat(batch)
  analysed <- list(whole_grid = strat_error(estimates),
                   survey_frame = strat_error(estimates, scope = "survey_frame"))
  expect_equal(batch, original)
  expect_equal(sort(unique(estimates$total_strat$survey)), c("east", "west"))

  key <- population$sim_length(population$ages, length_age_key = TRUE)
  weight <- population$sim_weight(population$lengths)
  for (label in names(batch$survey_programs)) {
    view <- batch |> select_survey(label) |> run_strat()
    for (scope in names(analysed)) {
      single <- strat_error(view, scope = scope)
      restored <- select_survey(analysed[[scope]], label)
      expect_identical(data.table::key(analysed[[scope]]$age_strat_error),
                       c("survey", "sim", "year", "age"))
      for (field in c("total_strat", "length_strat", "age_strat", "survey_coverage",
                      "total_strat_error", "length_strat_error", "age_strat_error",
                      "total_strat_error_stats", "length_strat_error_stats", "age_strat_error_stats",
                      "length_strat_weight_error_stats", "age_strat_weight_error_stats")) {
        expect_equal(restored[[field]], single[[field]])
      }
      for (rep in 1:2) {
        for (yr in population$years) {
          frame <- single$survey_frame
          cells <- frame$cell[frame$sim == rep & frame$year == yr]
          spatial <- single$sp_N[single$sp_N$cell %in% cells & single$sp_N$year == yr, ]
          frame_age <- vapply(population$ages, function(a) sum(spatial$N[spatial$age == a]),
                              numeric(1))
          stock_age <- population$N[, as.character(yr)]
          joint <- sweep(key, 2, if (scope == "whole_grid") stock_age else frame_age, "*")
          catchable <- if (label == "west") sweep(joint, 2, age_q(population$ages), "*") else
            joint * length_q(population$lengths)
          total <- single$total_strat_error[sim == rep & year == yr]
          ages <- single$age_strat_error[sim == rep & year == yr & !is.na(age)]
          lengths <- single$length_strat_error[sim == rep & year == yr]
          expect_equal(total$N_frame, sum(frame_age))
          expect_equal(total$N_stock, sum(stock_age))
          expect_equal(total$I, sum(catchable))
          expect_equal(ages$I, unname(colSums(catchable)))
          expect_equal(lengths$I[match(population$lengths, lengths$length)],
                       unname(rowSums(catchable)))
          available_weight_at_age <- colSums(catchable * weight) / colSums(catchable)
          available_weight_at_age[is.nan(available_weight_at_age)] <- NA_real_
          expect_equal(ages$mean_weight_true_kg, unname(available_weight_at_age))
          expect_equal(lengths$mean_weight_true_kg[match(population$lengths, lengths$length)],
                       ifelse(rowSums(catchable) > 0, weight, NA_real_) |> unname())
          if (!length(cells)) {
            expect_true(is.na(total$I_hat))
            expect_identical(total$estimate_status, "no_sets")
          }
        }
      }
    }
  }
  expect_null(run_strat(analysed$survey_frame)$age_strat_error)
  expect_null(run_strat(analysed$survey_frame)$age_strat_weight_error_stats)
})

test_that("the default uses stored whole-grid truth when the frame is unknown", {
  withr::local_seed(1251)
  population <- small_population()
  sets <- sim_sets(population, x <= 0, min_sets = 1, set_den = 0, trawl_dim = c(1, 5))
  attr(sets, "survey_frame") <- NULL
  legacy <- population |>
    sim_survey(custom_sets = sets, q_by = "length",
                q = function(length) ifelse(length < 30, 0.8, 0.25)) |>
    run_strat()
  legacy[c("sp_N", "survey_selection", "survey_replicates")] <- NULL
  whole <- strat_error(legacy)
  framed <- strat_error(legacy, scope = "survey_frame")
  expect_equal(whole, strat_error(legacy, scope = "whole_grid"))
  expect_equal(whole$total_strat_error$I, unname(colSums(legacy$I)))
  expect_equal(whole$total_strat_error$error,
               whole$total_strat_error$I_hat - whole$total_strat_error$I)
  expect_true(all(whole$total_strat_error$truth_status == "ok"))
  expect_true(all(is.na(whole$total_strat_error$N_frame)))
  expect_true(all(is.na(whole$age_strat_error$mean_weight_true_kg)))
  expect_true(all(is.na(whole$length_strat_weight_error_stats)))
  expect_equal(whole$total_strat_error$N_stock, unname(colSums(legacy$N)))
  for (yr in legacy$years) {
    ages <- whole$age_strat_error[year == yr & !is.na(age)]
    lengths <- whole$length_strat_error[year == yr & length %in% legacy$lengths]
    expect_equal(ages$I, unname(legacy$I[as.character(ages$age), as.character(yr)]))
    expect_equal(lengths$I, unname(legacy$I_at_length[as.character(lengths$length), as.character(yr)]))
  }
  expect_equal(framed$total_strat_error$I_hat, whole$total_strat_error$I_hat)
  expect_true(all(is.finite(framed$total_strat_error$I_hat)))
  expect_true(all(is.na(framed$total_strat_error$I)))
  expect_true(all(is.na(framed$total_strat_error$error)))
  expect_true(all(framed$total_strat_error$truth_status == "frame_unknown"))
  expect_equal(framed$total_strat_error$N_stock, whole$total_strat_error$N_stock)
  expect_true(all(is.na(framed$total_strat_error_stats)))
  expect_error(strat_error(legacy, scope = "visited_cells"), "whole_grid.*survey_frame")

  legacy[["I"]] <- NULL
  length_only <- strat_error(legacy)
  expect_true(all(length_only$total_strat_error$truth_status == "truth_unavailable"))
  expect_true(all(length_only$age_strat_error[!is.na(age)]$truth_status == "truth_unavailable"))
  expect_equal(length_only$length_strat_error, whole$length_strat_error)
})

test_that("both scopes compare lengths only at compatible bin widths", {
  withr::local_seed(1252)
  survey <- small_survey()
  expect_warning(incompatible <- run_strat(survey, length_group = 1), "Matching length truth")
  coarse <- run_strat(survey, length_group = 6)
  for (scope in c("whole_grid", "survey_frame")) {
    mismatch <- strat_error(incompatible, scope = scope)
    expect_true(all(is.na(mismatch$length_strat_error$I)))
    expect_true(all(is.na(mismatch$length_strat_error$mean_weight_true_kg)))
    expect_true(all(mismatch$length_strat_error$truth_status == "length_bins_mismatch"))
    expect_true(all(is.finite(mismatch$total_strat_error$I)))
    compared <- strat_error(coarse, scope = scope)
    sums <- compared$length_strat_error[, list(I_length = sum(I, na.rm = TRUE)), by = "year"]
    expect_equal(sums$I_length, compared$total_strat_error$I)
    native <- data.table::data.table(length = group_lengths(survey$lengths, 6),
                                    I = survey$I_at_length[, 1],
                                    kg = survey$sim_weight(survey$lengths))
    expected <- native[I > 0, list(weight = stats::weighted.mean(kg, I)), by = "length"]
    result <- compared$length_strat_error[year == survey$years[1]]
    expect_equal(result$mean_weight_true_kg[match(expected$length, result$length)], expected$weight)
    finite <- compared$length_strat_error[is.finite(mean_weight_error_kg)]
    expect_equal(finite$mean_weight_error_kg, finite$mean_weight_kg - finite$mean_weight_true_kg)
    expect_equal(compared$length_strat_weight_error_stats, error_stats(finite$mean_weight_error_kg))
  }
  incompatible$survey_frame <- NULL
  whole <- strat_error(incompatible)
  framed <- strat_error(incompatible, scope = "survey_frame")
  expect_true(all(whole$length_strat_error$truth_status == "length_bins_mismatch"))
  expect_true(all(framed$length_strat_error$truth_status == "frame_unknown"))
})

test_that("zeros, missing composition and partial age allocation are distinct", {
  withr::local_seed(1321)
  zero <- small_population() |>
    sim_survey(q = function(x) 0, set_den = 0, min_sets = 1) |>
    run_strat() |>
    strat_error()
  expect_true(all(zero$total_strat$total == 0 & zero$total_strat$estimate_status == "ok"))
  expect_true(all(zero$length_strat$total == 0))
  expect_true(all(zero$age_strat[!is.na(age)]$total == 0))

  no_lengths <- small_survey(lengths_cap = 0, ages_cap = 0, binom_error = FALSE) |>
    run_strat() |>
    strat_error()
  expect_true(all(no_lengths$total_strat$total > 0))
  expect_true(all(is.na(no_lengths$length_strat$total)))
  expect_true(all(is.na(no_lengths$age_strat$total)))
  expect_true(all(no_lengths$length_strat$estimate_status == "missing_lengths"))

  partial <- small_survey()
  partial$setdet$n <- rep(4, nrow(partial$setdet))
  partial$samp <- data.table::data.table(set = rep(partial$setdet$set, each = 2),
                                         age = c(1, 99), length = c(10, 20),
                                         measured = TRUE, aged = c(TRUE, FALSE))
  partial <- partial |> run_strat() |> strat_error()
  expect_equal(partial$age_strat[is.na(age)]$total, partial$total_strat$total / 2)
  expect_equal(partial$age_strat[age == 1]$total, partial$total_strat$total / 2)
  expect_true(all(partial$age_strat[!is.na(age)]$estimate_status == "partial_age"))
  expect_true(all(is.na(partial$age_strat[!is.na(age) & age != 1]$total)))
  expect_true(all(is.na(partial$age_strat_error$error)))
  expect_true(all(is.na(partial$age_strat_error[is.na(age)]$I)))
  expect_true(all(is.na(partial$age_strat_error_stats)))
  partial$samp$aged[] <- FALSE
  unaged <- run_strat(partial)
  expect_equal(unaged$age_strat[is.na(age)]$total, unaged$total_strat$total)
  expect_true(all(is.na(unaged$age_strat[!is.na(age)]$total)))
  expect_true(all(is.finite(unaged$length_strat$total)))
})

test_that("missing intended strata suppress full-frame estimates", {
  withr::local_seed(1374)
  grid <- make_grid(x_range = c(-10, 10), y_range = c(-10, 10), res = c(5, 5),
                      shelf_width = 4, shelf_depth = 50, depth_range = c(0, 100),
                      strat_breaks = c(0, 101), strat_splits = 2, method = "linear")
  population <- sim_distribution(small_abundance(), grid = grid, ays_covar = zero_covariance)
  sets <- sim_sets(population, set_den = 0, min_sets = 1, trawl_dim = c(1, 5))
  frame <- attr(sets, "survey_frame")
  omitted <- sets$strat[1]
  sets <- sets[sets$strat != omitted]
  expect_gt(nrow(sets), 0)
  result <- population |>
    sim_survey(custom_sets = sets, survey_frame = frame, q = function(x) 1) |>
    run_strat() |>
    strat_error(scope = "survey_frame")
  expect_true(all(result$survey_coverage[strat == omitted]$sets_sampled == 0))
  expect_true(all(result$total_strat$estimate_status == "unsampled_strata"))
  expect_true(all(is.na(result$total_strat$total)))
  expect_true(all(is.na(result$length_strat$total)))
  expect_true(all(is.na(result$age_strat$total)))
  expect_true(all(result$length_strat$weight_status == "unsampled_strata"))
  expect_true(all(result$age_strat[!is.na(age)]$maturity_status == "unsampled_strata"))
  expect_equal(result$total_strat_error$I, unname(colSums(result$N)))
})

test_that("survey estimates use observations without spatial or latent catch truth", {
  withr::local_seed(1444)
  population <- small_population()
  config <- sim_survey(lengths_cap = 3, ages_cap = 1, set_den = 0,
                         min_sets = 1, trawl_dim = c(1, 5), light = FALSE)
  batch <- sim_surveys(population, first = config, second = config)
  estimates <- run_strat(batch)
  for (label in names(batch$survey_programs)) {
    program <- batch$survey_programs[[label]]
    program$samp <- data.table::copy(program$samp[program$samp$measured, ])
    program$samp$age[!program$samp$aged] <- 99
    program[c("sp_I", "I", "I_at_length", "sample_truth", "full_setdet")] <- NULL
    batch$survey_programs[[label]] <- program
  }
  batch[c("sp_N", "N", "N_at_length")] <- NULL
  observed_only <- run_strat(batch)
  expect_equal(observed_only$total_strat, estimates$total_strat)
  expect_equal(observed_only$length_strat, estimates$length_strat)
  expect_equal(observed_only$age_strat, estimates$age_strat)
})

test_that("scenario fan intervals tolerate missing estimates", {
  comparisons <- data.table::data.table(year = c(2001, 2001, 2001, 2002),
                                         I_hat = c(1, NA, 3, NA))
  intervals <- .ints(comparisons, quants = 50, by = "year")
  expect_equal(intervals$lower, c(1.5, NA))
  expect_equal(intervals$upper, c(2.5, NA))
})
