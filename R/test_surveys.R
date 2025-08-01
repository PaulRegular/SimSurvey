

#' Set up a series of surveys from all combinations of settings supplied
#'
#' A convenience function that wraps [`base::expand.grid()`] to generate all combinations
#' of survey design parameters and adds a unique survey number to each row.
#'
#' @param set_den A vector of set densities (sets per grid unit squared).
#' @param lengths_cap A vector of maximum numbers of lengths measured per set.
#' @param ages_cap A vector of maximum numbers of otoliths to collect per length group
#' per division per year.
#'
#' @return A `data.frame` containing all combinations of the supplied vectors, with an added
#' `survey` column identifying each combination.
#'
#' @export

expand_surveys <- function(set_den = c(0.5, 1, 2, 5, 10) / 1000,
                           lengths_cap = c(5, 10, 20, 50, 100, 500, 1000),
                           ages_cap = c(2, 5, 10, 20, 50)) {
  surveys <- expand.grid(set_den = set_den, lengths_cap = lengths_cap, ages_cap = ages_cap)
  data.table(survey = seq.int(nrow(surveys)), surveys)
}



#' Test sampling design of multiple surveys using a stratified analysis
#'
#' This function allows a series of sampling design settings to be tested on a
#' simulated population. True population values are compared to stratified estimates
#' of abundance using a user-specified number of simulated surveys.
#'
#' @param sim A simulation object returned by [`sim_distribution()`].
#' @param surveys A `data.frame` or `data.table` of survey configurations, formatted like
#' the object returned by [`expand_surveys()`].
#' @param keep_details Integer. Retain full details for one survey (specified by survey number),
#' and drop the rest to reduce object size.
#' @param n_sims Number of surveys to simulate per design. Large values may consume significant RAM.
#' @param n_loops Number of times to loop [`sim_survey()`].
#' Total number of simulations = `n_sims` × `n_loops`.
#' A lower `n_sims` and higher `n_loops` combination is more memory efficient but may take longer.
#' @param cores Number of processor cores to use in parallel.
#' @param export_dir Optional directory path to export intermediate results.
#' Useful for resuming later with [`resume_test()`]. If `NULL`, nothing is exported.
#' @param progress Logical. Should progress bar and messages be displayed?
#' @inheritParams run_strat
#' @inheritDotParams sim_survey -sim -n_sims -set_den -lengths_cap -ages_cap -light
#'
#' @details
#' Depending on the number of surveys and simulations, `test_surveys()` can take a long time to run.
#'
#' The [`resume_test()`] function can be used to resume partial runs.
#' Note: progress bar time estimates may be biased if resuming previously completed iterations.
#'
#' Internally, this function calls a helper called `test_loop()` to process each survey simulation.
#'
#' **Caution:** When using `...` inside [`resume_test()`], be careful not to pass arguments that
#' were not part of the original `test_surveys()` call, as this could change simulation settings.
#'
#' @return The returned object includes:
#'
#' - A table of survey designs tested
#' - Stratified error results (`*_strat_error` and `*_strat_error_stats`)
#' - Error statistics:
#'   - `ME`: Mean error
#'   - `MAE`: Mean absolute error
#'   - `MSE`: Mean squared error
#'   - `RMSE`: Root mean squared error
#' - A summary table of total sample sizes (`samp_totals`)
#'
#' Survey and stratified analysis details are dropped for all but one retained survey (via `keep_details`).
#'
#' @examples
#' \donttest{
#' pop <- sim_abundance(ages = 1:20, years = 1:5) %>%
#'   sim_distribution(grid = make_grid(res = c(10, 10)))
#'
#' surveys <- expand_surveys(
#'   set_den = c(1, 2) / 1000,
#'   lengths_cap = c(100, 500),
#'   ages_cap = c(5, 20)
#' )
#'
#' # Simulate 25 surveys for each of 8 survey designs (low for example speed)
#' tests <- test_surveys(
#'   pop, surveys = surveys, keep_details = 1,
#'   n_sims = 5, n_loops = 5, cores = 1
#' )
#'
#' library(plotly)
#' tests$total_strat_error %>%
#'   filter(survey == 8, sim %in% 1:50) %>%
#'   group_by(sim) %>%
#'   plot_ly(x = ~year) %>%
#'   add_lines(y = ~I_hat, alpha = 0.5, name = "estimated") %>%
#'   add_lines(y = ~I, color = I("black"), name = "true") %>%
#'   layout(xaxis = list(title = "Year"),
#'          yaxis = list(title = "Abundance index"))
#'
#' plot_total_strat_fan(tests, surveys = 1:8)
#' plot_length_strat_fan(tests, surveys = 1:8)
#' plot_age_strat_fan(tests, surveys = 1:8)
#' plot_age_strat_fan(tests, surveys = 1:8, select_by = "age")
#'
#' plot_error_surface(tests, plot_by = "rule")
#' plot_error_surface(tests, plot_by = "samples")
#'
#' plot_survey_rank(tests, which_strat = "length")
#' plot_survey_rank(tests, which_strat = "age")
#' }
#'
#' @export
#'
#' @import progress
#' @import doParallel
#' @import parallel
#' @import foreach

test_surveys <- function(sim, surveys = expand_surveys(), keep_details = 1,
                         n_sims = 1, n_loops = 100, cores = 2, export_dir = NULL,
                         length_group = "inherit", alk_scale = "division",
                         progress = TRUE, ...) {

  if (!is.null(export_dir)) {
    save(list = ls(all.names = TRUE),
         file = file.path(export_dir, "test_inputs.RData"))
  }
  .test_loop(sim = sim, surveys = surveys, n_sims = n_sims,
             n_loops = n_loops, cores = cores, export_dir = export_dir,
             complete = NULL, keep_details = keep_details,
             length_group = length_group, alk_scale = alk_scale,
             progress = progress, ...)

}


#' @export
#' @rdname test_surveys
resume_test <- function(export_dir = NULL, ...) {
  sim <- surveys <- n_sims <- n_loops <- cores <-
    complete <- keep_details <- length_group <- alk_scale <-
    progress <- NULL
  load(file.path(export_dir, "test_inputs.RData"))
  load(file.path(export_dir, "complete.RData"))
  .test_loop(sim = sim, surveys = surveys, n_sims = n_sims,
             n_loops = n_loops, cores = cores, export_dir = export_dir,
             complete = complete, keep_details = keep_details,
             length_group = length_group, alk_scale = alk_scale,
             progress = progress, ...)
}



## Helper function for test_surveys and resume_test
.test_loop <- function(sim = NULL, surveys = NULL, n_sims = NULL,
                       n_loops = NULL, cores = NULL, export_dir = NULL,
                       complete = NULL, keep_details = NULL,
                       length_group = NULL, alk_scale = NULL, progress = NULL,
                       ...) {

  j <- error <- NULL

  ## Containers
  samp_totals <- total_strat_error <- length_strat_error <- age_strat_error <- vector("list", nrow(surveys))
  if (is.null(complete) && !is.null(export_dir)) {
    dir.create(file.path(export_dir, "samp_totals"), showWarnings = FALSE)
    dir.create(file.path(export_dir, "total_strat_error"), showWarnings = FALSE)
    dir.create(file.path(export_dir, "length_strat_error"), showWarnings = FALSE)
    dir.create(file.path(export_dir, "age_strat_error"), showWarnings = FALSE)
  }
  if (is.null(complete)) {
    complete <- rep(FALSE, nrow(surveys))
  }

  if (any(!complete)) {

    incomplete <- which(!complete)

    if (progress) {
      message("\nRunning simulations...\n")
      pb <- progress::progress_bar$new(
        format = "Survey :current of :total [:bar] :percent in :elapsed (eta: :eta)",
        total = nrow(surveys), clear = FALSE, show_after = 0, width = 100)
      invisible(pb$tick(min(incomplete) - 1))
    }

    for (i in incomplete) {

      cl <- parallel::makeCluster(cores) # use parallel computation
      doParallel::registerDoParallel(cl, cores = cores)
      loop_error <- foreach::foreach(j = seq(n_loops), .packages = "SimSurvey") %dopar% {
                              res <- sim_survey(sim, n_sims = n_sims,
                                                set_den = surveys$set_den[i],
                                                lengths_cap = surveys$lengths_cap[i],
                                                ages_cap = surveys$ages_cap[i],
                                                ...) %>%
                                run_strat(length_group = length_group,
                                          alk_scale = alk_scale) %>%
                                strat_error()
                              samp_totals <- res$samp_totals
                              samp_totals$sim <- samp_totals$sim + (j * n_sims - n_sims)
                              total_strat_error <- res$total_strat_error
                              total_strat_error$sim <- total_strat_error$sim + (j * n_sims - n_sims)
                              length_strat_error <- res$length_strat_error
                              length_strat_error$sim <- length_strat_error$sim + (j * n_sims - n_sims)
                              age_strat_error <- res$age_strat_error
                              age_strat_error$sim <- age_strat_error$sim + (j * n_sims - n_sims)
                              list(samp_totals = samp_totals,
                                   total_strat_error = total_strat_error,
                                   length_strat_error = length_strat_error,
                                   age_strat_error = age_strat_error)
                            }
      parallel::stopCluster(cl) # stop parallel process

      samp_totals_i <- data.table::rbindlist(lapply(loop_error, `[[`, "samp_totals"))
      total_strat_error_i <- data.table::rbindlist(lapply(loop_error, `[[`, "total_strat_error"))
      length_strat_error_i <- data.table::rbindlist(lapply(loop_error, `[[`, "length_strat_error"))
      age_strat_error_i <- data.table::rbindlist(lapply(loop_error, `[[`, "age_strat_error"))
      samp_totals_i$survey <- i
      total_strat_error_i$survey <- i
      length_strat_error_i$survey <- i
      age_strat_error_i$survey <- i

      if (is.null(export_dir)) {
        samp_totals[[i]] <- samp_totals_i
        total_strat_error[[i]] <- total_strat_error_i
        length_strat_error[[i]] <- length_strat_error_i
        age_strat_error[[i]] <- age_strat_error_i
      } else {
        data.table::fwrite(samp_totals_i, file = file.path(export_dir, "samp_totals", paste0(i, ".csv")))
        data.table::fwrite(total_strat_error_i, file = file.path(export_dir, "total_strat_error", paste0(i, ".csv")))
        data.table::fwrite(length_strat_error_i, file = file.path(export_dir, "length_strat_error", paste0(i, ".csv")))
        data.table::fwrite(age_strat_error_i, file = file.path(export_dir, "age_strat_error", paste0(i, ".csv")))
        complete[i] <- TRUE
        save(complete, file = file.path(export_dir, "complete.RData"))
      }

      if (progress) pb$tick()

    }

  }

  ## Compile results
  if (progress) message("\nCompiling results...\n")
  if (!is.null(export_dir)) {
    if (progress) {
      pb <- progress::progress_bar$new(
        format = "Survey :current of :total [:bar] :percent in :elapsed (eta: :eta)",
        total = nrow(surveys), clear = FALSE, show_after = 0, width = 100)
    }
    for (i in seq(nrow(surveys))) {
      samp_totals[[i]] <- data.table::fread(file.path(export_dir, "samp_totals", paste0(i, ".csv")),
                                            verbose = FALSE, showProgress = FALSE)
      total_strat_error[[i]] <- data.table::fread(file.path(export_dir, "total_strat_error", paste0(i, ".csv")),
                                                  verbose = FALSE, showProgress = FALSE)
      length_strat_error[[i]] <- data.table::fread(file.path(export_dir, "length_strat_error", paste0(i, ".csv")),
                                                   verbose = FALSE, showProgress = FALSE)
      age_strat_error[[i]] <- data.table::fread(file.path(export_dir, "age_strat_error", paste0(i, ".csv")),
                                                verbose = FALSE, showProgress = FALSE)
      if (progress) pb$tick()
    }
  }
  sim$surveys <- data.table::data.table(surveys)
  sim$samp_totals <- data.table::rbindlist(samp_totals); rm(samp_totals)
  sim$total_strat_error <- data.table::rbindlist(total_strat_error); rm(total_strat_error)
  sim$length_strat_error <- data.table::rbindlist(length_strat_error); rm(length_strat_error)
  sim$age_strat_error <- data.table::rbindlist(age_strat_error); rm(age_strat_error)

  ## Calculate some stats
  sim$total_strat_error_stats <- sim$total_strat_error[, as.list(error_stats(error)), by = "survey"]
  sim$length_strat_error_stats <- sim$length_strat_error[, as.list(error_stats(error)), by = "survey"]
  sim$age_strat_error_stats <- sim$age_strat_error[, as.list(error_stats(error)), by = "survey"]

  ## Keep details from one survey
  i <- which(surveys$survey == keep_details)
  res <- sim_survey(sim, n_sims = n_sims,
                    set_den = surveys$set_den[i],
                    lengths_cap = surveys$lengths_cap[i],
                    ages_cap = surveys$ages_cap[i],
                    light = FALSE,
                    ...) %>% run_strat()
  sim$N0 <- res$N[, 1] # attach rounded results
  sim$R <- res$N[1, ]
  sim$N <- res$N
  sim$sp_N <- res$sp_N
  sim$I <- res$I
  sim$I_at_length <- res$length$I
  sim$full_setdet <- res$full_setdet
  sim$setdet <- res$setdet
  sim$samp <- res$samp
  sim$total_strat <- res$total_strat
  sim$length_strat <- res$length_strat
  sim$age_strat <- res$age_strat

  if (!is.null(export_dir)) {
    save(sim, file = file.path(export_dir, "test_output.RData"))
  }

  sim

}




