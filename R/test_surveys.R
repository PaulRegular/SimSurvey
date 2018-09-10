

#' Set-up a series of surveys from all combinations of settings supplied
#'
#' @description Function is simply a wrapper for \code{\link{base::expand.grid}} that
#' adds a survey number to the returned object
#'
#' @param set_den           Vector of set densities (number of sets per [grid unit] squared)
#' @param lengths_cap       Vector of maximum number of lengths measured per set
#' @param ages_cap          Vector of maximum number of otoliths to collect per length group
#'                          per division per year
#'
#' @export
#'

expand_surveys <- function(set_den = c(0.5, 1, 2, 5, 10) / 1000,
                           lengths_cap = c(5, 10, 20, 50, 100, 500, 1000),
                           ages_cap = c(2, 5, 10, 20, 50)) {
  surveys <- expand.grid(set_den = set_den, lengths_cap = lengths_cap, ages_cap = ages_cap)
  data.table(survey = seq.int(nrow(surveys)), surveys)
}




## helper functions for test_surveys and resume_test
.test_loop <- function(sim = NULL, surveys = NULL, n_sims = NULL,
                       n_loops = NULL, cores = NULL, export = NULL,
                       export_dir = NULL, survey_error = NULL, ...) {

  if (is.null(survey_error)) {
    survey_error <- vector("list", nrow(surveys))
  }
  incomplete <- which(sapply(survey_error, is.null))

  pb <- progress::progress_bar$new(
    format = "Survey :current of :total [:bar] :percent in :elapsed (eta: :eta)",
    total = nrow(surveys), clear = FALSE, show_after = 0, width = 100)
  invisible(pb$tick(min(incomplete) - 1))

  for (i in incomplete) {

    cl <- makeCluster(cores) # use parallel computation
    registerDoParallel(cl)
    loop_error <- foreach(j = seq(n_loops),
                          .packages = "SimSurvey") %dopar% {
                            res <- sim_survey(sim, n_sims = n_sims,
                                              set_den = surveys$set_den[i],
                                              lengths_cap = surveys$lengths_cap[i],
                                              ages_cap = surveys$ages_cap[i],
                                              ...) %>%
                              run_strat() %>% strat_error()
                            samp_totals <- res$samp_totals
                            samp_totals$sim <- samp_totals$sim + (j * n_sims - n_sims)
                            total_strat_error <- res$total_strat_error
                            total_strat_error$sim <- total_strat_error$sim + (j * n_sims - n_sims)
                            age_strat_error <- res$age_strat_error
                            age_strat_error$sim <- age_strat_error$sim + (j * n_sims - n_sims)
                            list(samp_totals = samp_totals,
                                 total_strat_error = total_strat_error,
                                 age_strat_error = age_strat_error)
                          }
    stopCluster(cl) # stop parallel process

    samp_totals <- data.table::rbindlist(lapply(loop_error, `[[`, "samp_totals"))
    total_strat_error <- data.table::rbindlist(lapply(loop_error, `[[`, "total_strat_error"))
    age_strat_error <- data.table::rbindlist(lapply(loop_error, `[[`, "age_strat_error"))
    samp_totals$survey <- i
    total_strat_error$survey <- i
    age_strat_error$survey <- i
    survey_error[[i]]$samp_totals <- samp_totals
    survey_error[[i]]$total_strat_error <- total_strat_error
    survey_error[[i]]$age_strat_error <- age_strat_error

    if (!is.null(export)) {
      save(survey_error, file = file.path(export_dir, "survey_error.RData"))
    }

    pb$tick()

  }

  survey_error

}

.test_fin <- function(sim = NULL,  n_sims = NULL, surveys = NULL, keep_details = NULL,
                      survey_error = NULL, export = NULL, export_dir = NULL, ...) {

  sim$surveys <- data.table::data.table(surveys)
  sim$samp_totals <- do.call(rbind, lapply(survey_error, `[[`, "samp_totals"))
  sim$total_strat_error <- data.table::rbindlist(lapply(survey_error, `[[`, "total_strat_error"))
  sim$age_strat_error <- data.table::rbindlist(lapply(survey_error, `[[`, "age_strat_error"))
  sim$total_strat_error_stats <- sim$total_strat_error[, list(MAE = mean(abs(error)),
                                                              MSE = mean(error ^ 2),
                                                              RMSE = sqrt(mean(error ^ 2))),
                                                       by = "survey"]
  sim$age_strat_error_stats <- sim$age_strat_error[, list(MAE = mean(abs(error)),
                                                          MSE = mean(error ^ 2),
                                                          RMSE = sqrt(mean(error ^ 2))),
                                                   by = "survey"]

  i <- which(surveys$survey == keep_details)
  res <- sim_survey(sim, n_sims = n_sims,
                    set_den = surveys$set_den[i],
                    lengths_cap = surveys$lengths_cap[i],
                    ages_cap = surveys$ages_cap[i],
                    light = FALSE,
                    ...) %>% run_strat()
  sim$I <- res$I
  sim$sp_I <- res$sp_I
  sim$full_setdet <- res$full_setdet
  sim$setdet <- res$setdet
  sim$samp <- res$samp
  sim$total_strat <- res$total_strat
  sim$length_strat <- res$length_strat
  sim$age_strat <- res$age_strat

  if (!is.null(export)) {
    save(sim, file = file.path(export_dir, "test_output.RData"))
  }

  sim

}



#' Test sampling design of multiple surveys
#'
#' @description   This function allows a series of sampling design settings to
#' be set and tested on the simulated population.
#'
#' @param sim               Simulation from \code{\link{sim_distribution}}.
#' @param surveys           A data.frame or data.table with a sequence of surveys and their settings
#'                          with a format like the data.table returned by \code{\link{expand_surveys}}.
#' @param keep_details      Keep details of a specific survey in the data.frame supplied to \code{surveys}.
#'                          Survey and stratified analysis details are dropped for the other surveys to
#'                          minimize object size.
#' @param n_sims            Number of times to simulate a survey over the simulated population.
#'                          Requesting a large number of simulations here may max out your RAM.
#' @param n_loops           Number of times to run the \code{\link{sim_survey}} function. Total
#'                          simulations run will be the product of \code{n_sims} and \code{n_loops}
#'                          arguments. Low numbers of \code{n_sims} and high numbers of \code{n_loops}
#'                          will be easier on RAM, but may be slower.
#' @param cores             Number of cores to use in parallel. More cores should speed up the process.
#' @param export            Directory for exporting results as they are generated. If NULL, nothing is
#'                          exported.
#' @param ...               Additional arguments to pass to \code{\link{sim_survey}}.
#'
#' @return Adds a table of survey designs tested. Also adds details and summary
#'         stats of stratified estimate error to the \code{sim} list, ending with
#'         \code{"_strat_error"} or \code{"_strat_error_stats"}. Error statistics
#'         includes mean absolute error (\code{"MAE"}), mean squared error (\code{"MSE"}),
#'         and root mean squared error (\code{"RMSE"}). Also adds a sample size summary table
#'         (\code{"samp_totals"}) to the list. Survey and stratified analysis
#'         details are not kept to minimize object size.
#'
#' @export
#'
#' @import progress
#' @import doParallel
#' @import parallel
#' @import foreach
#'

test_surveys <- function(sim, surveys = expand_surveys(), keep_details = 1,
                         n_sims = 1, n_loops = 100, cores = 2, export = NULL, ...) {

  if (!is.null(export)) {
    export_dir <- file.path(export, paste0(Sys.Date(), "_test"))
    dir.create(export_dir, showWarnings = FALSE)
    save(list = ls(all.names = TRUE),
         file = file.path(export_dir, "test_inputs.RData"))
  }

  survey_error <- .test_loop(sim = sim, surveys = surveys, n_sims = n_sims,
             n_loops = n_loops, cores = cores, export = export,
             export_dir = export_dir, survey_error = NULL, ...)
  .test_fin(sim = sim, surveys = surveys, survey_error = survey_error,
            n_sims = n_sims, keep_details = keep_details,
            export = export, export_dir = export_dir, ...)


}


#' Resume partial run of \code{\link{test_surveys}}
#'
#' @param dir  Export directory specified when \code{\link{test_surveys}} was run.
#'
#' @details Progress bar time estimates will be biased here by previous completions
#'
#' @export
#'

resume_test <- function(dir = NULL) {
  load(file.path(dir, "test_inputs.RData"))
  load(file.path(dir, "survey_error.RData"))
  survey_error <- .test_loop(sim = sim, surveys = surveys, n_sims = n_sims,
                             n_loops = n_loops, cores = cores, export = export,
                             export_dir = export_dir, survey_error = survey_error, ...)
  .test_fin(sim = sim, surveys = surveys, survey_error = survey_error,
            n_sims = n_sims, keep_details = keep_details,
            export = export, export_dir = export_dir, ...)
}


