test_that("sets retain eligible cells before allocation by replicate and year", {
  withr::local_seed(861)
  pop <- small_population()
  sets <- sim_sets(pop, year != 2002 & cell <= ifelse(sim == 1, 8, 4),
                     n_sims = 2, min_sets = 2, set_den = 0, trawl_dim = c(1, 5))
  frame <- attr(sets, "survey_frame")
  expected <- data.table::CJ(sim = 1:2, year = c(2001, 2003), cell = 1:8)
  expected <- expected[sim == 1 | cell <= 4]
  expect_equal(as.data.frame(frame), as.data.frame(expected))
  expect_lt(nrow(sets), nrow(frame))
  expect_error(sim_sets(pop, NA), "without NAs")

  result <- sim_survey(pop, n_sims = 2, custom_sets = sets, q = function(age) 0.5)
  expect_equal(result$survey_frame, frame)
  expect_identical(result$survey_replicates, 1:2)
  expect_equal(.survey_frame_strata(result)$eligible_cells, c(8L, 8L, 4L, 4L))
  expect_equal(result$survey_selection$values,
               stats::setNames(rep(0.5, 3), as.character(pop$ages)))
  expect_equal(result$I, result$N * 0.5)
})

test_that("custom sets declare a frame explicitly or retain an unknown frame", {
  withr::local_seed(862)
  pop <- small_population()
  sets <- sim_sets(pop, cell <= 8, min_sets = 2, set_den = 0, trawl_dim = c(1, 5))
  frame <- attr(sets, "survey_frame")
  bare <- data.table::copy(sets)
  data.table::setattr(bare, "survey_frame", NULL)

  set.seed(863)
  unknown <- sim_survey(pop, custom_sets = bare)
  rng <- .Random.seed
  set.seed(863)
  declared <- sim_survey(pop, custom_sets = bare, survey_frame = frame)
  expect_null(unknown$survey_frame)
  expect_equal(declared$survey_frame, frame)
  expect_equal(unknown$setdet, declared$setdet)
  expect_equal(unknown$samp, declared$samp)
  expect_identical(.Random.seed, rng)
  expect_null(attr(bare, "survey_frame"))

  expect_error(sim_survey(pop, survey_frame = frame), "with custom_sets")
  expect_error(sim_survey(pop, custom_sets = bare, survey_frame = rbind(frame, frame[1])),
               "unique, nonmissing")
  invalid_cell <- data.table::copy(frame)
  invalid_cell$cell[1] <- max(pop$grid_xy$cell) + 1
  expect_error(sim_survey(pop, custom_sets = bare, survey_frame = invalid_cell),
               "replicate, year, and cell identifiers")
  outside <- frame[!bare[1], on = c("sim", "year", "cell")]
  expect_error(sim_survey(pop, custom_sets = bare, survey_frame = outside),
               "must lie within survey_frame")
  wrong_area <- data.table::copy(bare)
  wrong_area$strat_area[1] <- 2 * wrong_area$strat_area[1]
  expect_error(sim_survey(pop, custom_sets = wrong_area, survey_frame = frame),
               "strat_area must match")
})

test_that("survey frames remain known when no sets are allocated", {
  withr::local_seed(864)
  pop <- small_population()
  no_sets <- sim_sets(pop, min_sets = 0, set_den = 0)
  expect_equal(nrow(no_sets), 0L)
  full <- sim_survey(pop, min_sets = 0, set_den = 0)
  expect_equal(full$survey_frame, attr(no_sets, "survey_frame"))
  expect_equal(nrow(full$survey_frame), nrow(pop$grid_xy) * length(pop$years))
  expect_equal(nrow(full$setdet), 0L)
  expect_equal(nrow(full$samp), 0L)
  empty_sets <- sim_sets(pop, FALSE)
  empty <- sim_survey(pop, custom_sets = empty_sets)
  expect_s3_class(empty$survey_frame, "data.table")
  expect_equal(nrow(empty$survey_frame), 0L)
  expect_identical(empty$survey_replicates, 1L)
})

test_that("named surveys carry independent frames and evaluated catchability", {
  withr::local_seed(865)
  pop <- small_population()
  west <- sim_sets(pop, cell <= 8, min_sets = 2, set_den = 0, trawl_dim = c(1, 5))
  east <- sim_sets(pop, cell >= 5, min_sets = 2, set_den = 0, trawl_dim = c(1, 5))
  batch <- sim_surveys(pop,
    west = sim_survey(custom_sets = west, q = function(age) age >= 2),
    east = sim_survey(custom_sets = east, q_by = "length",
                      q = function(length) as.numeric(length >= 30 & length < 60))
  )
  expect_null(batch$survey_frame)
  a <- select_survey(batch, "west")
  b <- select_survey(batch, "east")
  expect_equal(a$survey_frame, attr(west, "survey_frame"))
  expect_equal(b$survey_frame, attr(east, "survey_frame"))
  expect_identical(a$survey_selection$by, "age")
  expect_equal(unname(a$survey_selection$values), c(0, 1, 1))
  expect_identical(b$survey_selection$by, "length")
  expect_equal(names(b$survey_selection$values), rownames(b$survey_selection$key))
  expect_equal(colSums(b$I), colSums(b$I_at_length))
  expect_null(.survey_view(batch, "east", restore_spatial = FALSE)$sp_N$I)
})

test_that("parallel surveys renumber frame replicates with sets", {
  withr::local_seed(866)
  local_worker_library()
  pop <- small_population()
  sets <- sim_sets(pop, cell <= 8, min_sets = 2, set_den = 0, trawl_dim = c(1, 5))
  result <- sim_survey_parallel(pop, n_sims = 1, n_loops = 2, cores = 1,
                                  custom_sets = sets, quiet = TRUE)
  expect_identical(result$survey_replicates, 1:2)
  expected <- data.table::CJ(sim = 1:2, year = pop$years, cell = 1:8)
  expect_equal(as.data.frame(result$survey_frame), as.data.frame(expected))
  expect_equal(nrow(result$setdet[!result$survey_frame, on = c("sim", "year", "cell")]), 0L)
})
