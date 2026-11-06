# Run with Rscript after installing this version of SimSurvey.
library(SimSurvey)
library(data.table)
setDTthreads(1L)

scan_counts <- function(operation_ids, sample) {
  list(n_measured = vapply(operation_ids, function(id) {
    sum(sample$operation == id)
  }, numeric(1)),
  n_aged = vapply(operation_ids, function(id) {
    sum(sample$aged[sample$operation == id])
  }, numeric(1)),
  n_maturity = vapply(operation_ids, function(id) {
    sum(sample$maturity_observed[sample$operation == id])
  }, numeric(1)))
}

benchmark_counts <- function(n_operations, repeats = 3L, grouped_calls = 100L) {
  operation_ids <- seq_len(n_operations)
  # Every fifth operation contributes no samples.
  sample <- data.table(operation = rep(operation_ids[operation_ids %% 5L != 0L], each = 24L))
  sample$aged <- seq_len(nrow(sample)) %% 3L == 0L
  sample$maturity_observed <- sample$aged & seq_len(nrow(sample)) %% 2L == 0L

  grouped_counts <- SimSurvey:::.fishery_sample_counts
  stopifnot(identical(scan_counts(operation_ids, sample), grouped_counts(operation_ids, sample)))
  scan_seconds <- replicate(repeats, system.time(scan_counts(operation_ids, sample))[["elapsed"]])
  grouped_seconds <- replicate(repeats, system.time({
    for (i in seq_len(grouped_calls)) grouped_counts(operation_ids, sample)
  })[["elapsed"]] / grouped_calls)
  data.table(n_operations = n_operations, n_sample_rows = nrow(sample),
             scan_seconds = median(scan_seconds), grouped_seconds = median(grouped_seconds))
}

print(rbindlist(lapply(c(400L, 1600L, 4000L), benchmark_counts)))
