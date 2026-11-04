test_that("invalid age and year sequences fail before population simulation", {
  bad_sequences <- list(3:1, c(1, 1, 2), c(1, 3, 4), numeric(),
                         c(1, NA_real_), c(1, Inf), c("1", "2"))
  for (sequence in bad_sequences) {
    expect_error(sim_abundance(ages = sequence),
                 "Age and year sequences must be ascending and increment by 1")
    expect_error(sim_abundance(years = sequence),
                 "Age and year sequences must be ascending and increment by 1")
  }
})
