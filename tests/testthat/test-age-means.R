test_that("length-stratified means weight ages by expanded length counts", {
  totals <- data.frame(length = c(20, 40), total = c(90, 10))
  aged <- data.frame(age = factor(c(3, 3)), length = c(20, 40),
                     weight = c(1, 3), mature = c(FALSE, TRUE))
  out <- age_means(totals, aged, c("mature", "weight", "length"))
  expect_identical(key(out), c("age", "metric"))
  expect_equal(out$age, rep(3, 3))
  expect_equal(out$estimate, c(22, 0.1, 1.2))
  expect_identical(out$n, rep(2L, 3))
})

test_that("the paper's age-two maturity example gives one third", {
  totals <- data.frame(length = 1:5, total = c(360, 450, 540, 400, 70))
  aged_n <- c(36, 45, 54, 40, 7)
  age_two <- c(12, 13, 7, 1, 0)
  mature_two <- c(3, 4, 3, 1, 0)
  aged <- rbindlist(lapply(seq_along(aged_n), function(j) {
    data.table(length = j, age = rep(c(2, 3), c(age_two[j], aged_n[j] - age_two[j])),
               mature = seq_len(aged_n[j]) <= mature_two[j])
  }))
  expect_equal(age_means(totals, aged, "mature")[age == 2, estimate], 1 / 3)
})

test_that("strata pool by expansion and each metric has its own support", {
  totals <- data.frame(year = 1, division = c("a", "b", "b"), length = 30,
                       total = c(80, 8, 12))
  aged <- data.table(year = 1, division = rep(c("a", "b"), each = 2), length = 30,
                     age = rep(c(3, 5), 2), mature = c(0, NA, 1, 0), weight = c(1, 1, 3, 3))
  before <- copy(aged)
  out <- age_means(totals, aged, c("mature", "weight"), strata = c("year", "division"), by = "year")
  expect_identical(aged, before)
  expect_equal(out[age == 3 & metric == "mature", estimate], 0.2)
  expect_equal(out[age == 3 & metric == "weight", estimate], 1.4)
  expect_identical(out[age == 5 & metric == "mature", status], "missing_trait")
  expect_true(is.na(out[age == 5 & metric == "mature", estimate]))
  expect_equal(out[age == 5 & metric == "mature", n], 1)
})

test_that("unknown, unsupported and zero length counts are distinct", {
  totals <- data.frame(length = c(20, 40), total = c(90, 10))
  aged <- data.frame(age = 3, length = 20, mature = 0)
  out <- age_means(totals, aged, "mature", ages = 3:4)
  expect_identical(out$status, rep("no_aging_support", 2))
  expect_true(all(is.na(out$estimate)))
  expect_equal(out$n, c(1, 0))
  totals$total[2] <- 0
  out <- age_means(totals, aged, "mature", ages = 3:4)
  expect_equal(out$estimate, c(0, NA_real_))
  expect_identical(out$status, c("ok", "no_weighted_fish"))
  totals$total[2] <- NA_real_
  expect_identical(age_means(totals, aged, "mature")$status, "missing_length_total")
  expect_identical(age_means(totals[1, ], aged[0, ], "mature", ages = 3)$status,
                   "no_aging_support")
  expect_error(age_means(totals, transform(aged, length = 50), "mature"), "must match")
  expect_error(age_means(totals, transform(aged, mature = 2), "mature"), "Invalid")
  expect_error(age_means(totals, transform(aged, age = "adult"), "mature"), "Age must")
  expect_error(age_means(totals, aged, "mature", by = "year"), "subset")
})
