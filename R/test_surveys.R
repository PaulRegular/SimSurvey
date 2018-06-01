

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

expand_surveys <- function(set_den = c(0.3, 0.5, 0.8, 1, 2, 3, 6, 9) / 1000,
                           lengths_cap = c(2, 3, 5, 8, 10, 20, 30, 60, 90, 100, 200,
                                           400, 600, 1000),
                           ages_cap = c(2, 3, 5, 8, 10, 20, 30, 60)) {
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
                            total_strat_error <- res$total_strat_error
                            total_strat_error$sim <- total_strat_error$sim + (j * n_sims - n_sims)
                            age_strat_error <- res$age_strat_error
                            age_strat_error$sim <- age_strat_error$sim + (j * n_sims - n_sims)
                            list(total_strat_error = total_strat_error, age_strat_error = age_strat_error)
                          }
    stopCluster(cl) # stop parallel process

    total_strat_error <- data.table::rbindlist(lapply(loop_error, `[[`, "total_strat_error"))
    age_strat_error <- data.table::rbindlist(lapply(loop_error, `[[`, "age_strat_error"))
    total_strat_error$survey <- i
    age_strat_error$survey <- i
    survey_error[[i]]$total_strat_error <- total_strat_error
    survey_error[[i]]$age_strat_error <- age_strat_error
    survey_error[[i]]$total_strat_error_stats <- c(survey = i, error_stats(total_strat_error$error))
    survey_error[[i]]$age_strat_error_stats <- c(survey = i, error_stats(age_strat_error$error))

    if (!is.null(export)) {
      save(survey_error, file = file.path(export_dir, "survey_error.RData"))
    }

    pb$tick()

  }

  survey_error

}

.test_fin <- function(sim = NULL, surveys = NULL, survey_error = NULL) {

  total_strat_error <- data.table::rbindlist(lapply(survey_error, `[[`, "total_strat_error"))
  age_strat_error <- data.table::rbindlist(lapply(survey_error, `[[`, "age_strat_error"))
  total_strat_error_stats <- do.call(rbind, lapply(survey_error, `[[`, "total_strat_error_stats"))
  total_strat_error_stats <- data.table::data.table(total_strat_error_stats)
  age_strat_error_stats <- do.call(rbind, lapply(survey_error, `[[`, "age_strat_error_stats"))
  age_strat_error_stats <- data.table::data.table(age_strat_error_stats)

  sim$surveys <- data.table::data.table(surveys)
  sim$total_strat_error <- total_strat_error
  sim$total_strat_error_stats <- total_strat_error_stats
  sim$age_strat_error <- age_strat_error
  sim$age_strat_error_stats <- age_strat_error_stats
  sim

}



#' Test sampling design of multiple surveys
#'
#' @description   This function allows a series of sampling design settings to
#' be set and tested on the simulated population. All combinations of the
#' supplied settings (\code{set_den}, \code{lengths_cap}, \code{ages_cap}) are tested.
#'
#' @param sim               Simulation from \code{\link{sim_distribution}}.
#' @param surveys           A data.frame or data.table with a sequence of surveys and their settings
#'                          with a format like the data.table returned by \code{\link{expand_surveys}}.
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
#'         and root mean squared error (\code{"RMSE"}). Survey and stratified analysis
#'         details are not kept to minimize object size.
#'
#' @export
#'
#' @import progress
#' @import doParallel
#' @import parallel
#' @import foreach
#'

test_surveys <- function(sim, surveys = expand_surveys(), n_sims = 1,
                         n_loops = 100, cores = 2, export = NULL, ...) {

  if (!is.null(export)) {
    export_dir <- file.path(export, paste0(Sys.Date(), "_test"))
    dir.create(export_dir, showWarnings = FALSE)
    save(list = ls(all.names = TRUE),
         file = file.path(export_dir, "test_inputs.RData"))
  }

  survey_error <- .test_loop(sim = sim, surveys = surveys, n_sims = n_sims,
             n_loops = n_loops, cores = cores, export = export,
             export_dir = export_dir, survey_error = NULL, ...)
  .test_fin(sim = sim, surveys = surveys, survey_error = survey_error)


}


#' Resume partial run of \code{\link{test_surveys}}
#'
#' @param dir  Export directory specified when \code{\link{test_surveys}} was run.
#'
#' @export
#'

resume_test <- function(dir = NULL) {
  load(file.path(dir, "test_inputs.RData"))
  load(file.path(dir, "survey_error.RData"))
  survey_error <- .test_loop(sim = sim, surveys = surveys, n_sims = n_sims,
                             n_loops = n_loops, cores = cores, export = export,
                             export_dir = export_dir, survey_error = survey_error, ...)
  .test_fin(sim = sim, surveys = surveys, survey_error = survey_error)
}

