

#' Set-up a series of surveys from all combinations of settings supplied
#'
#' @description Function is simply a wrapper for \code{\link[base]{expand.grid}} that
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



#' Test sampling design of multiple surveys using a stratified analysis
#'
#' @description   This function allows a series of sampling design settings to
#' be set and tested on the simulated population. True population values are compared
#' to stratified estimates of abundance.
#'
#' @param sim               Simulation from \code{\link{sim_distribution}}.
#' @param surveys           A data.frame or data.table with a sequence of surveys and their settings
#'                          with a format like the data.table returned by \code{\link{expand_surveys}}.
#' @param keep_details      Survey and stratified analysis details are dropped here to minimize object
#'                          size. This argument allows the user to keep the details of one
#'                          survey by specifying the survey number in the data.frame supplied to \code{surveys}.
#' @param n_sims            Number of times to simulate a survey over the simulated population.
#'                          Requesting a large number of simulations here may max out your RAM.
#' @param n_loops           Number of times to run the \code{\link{sim_survey}} function. Total
#'                          simulations run will be the product of \code{n_sims} and \code{n_loops}
#'                          arguments. Low numbers of \code{n_sims} and high numbers of \code{n_loops}
#'                          will be easier on RAM, but may be slower.
#' @param cores             Number of cores to use in parallel. More cores should speed up the process.
#' @param export_dir        Directory for exporting results as they are generated. Main use of the export
#'                          is to allow this process to pick up where \code{test_survey} left off by
#'                          calling \code{resume_test}. If NULL, nothing is exported.
#' @inherit                 run_strat
#' @inheritDotParams        sim_survey
#'
#' @details Depening on the settings, \code{test_surveys} may take a long time to run.
#' The \code{resume_test} function is for resuming partial runs of \code{test_surveys}.
#' Note that progress bar time estimates will be biased here by previous completions.
#' \code{test_loop} is a helper function used in both \code{test_surveys} and
#' \code{resume_test}.
#'
#' @return Adds a table of survey designs tested. Also adds details and summary
#'         stats of stratified estimate error to the \code{sim} list, ending with
#'         \code{"_strat_error"} or \code{"_strat_error_stats"}. Error statistics
#'         includes mean absolute error (\code{"MAE"}), mean squared error (\code{"MSE"}),
#'         and root mean squared error (\code{"RMSE"}). Also adds a sample size summary table
#'         (\code{"samp_totals"}) to the list. Survey and stratified analysis
#'         details are not kept to minimize object size.
#'
#' @examples
#'
#' \dontrun{
#' pop <- sim_abundance(ages = 1:20, years = 1:5) %>%
#'            sim_distribution(grid = make_grid(res = c(10, 10)))
#'
#' surveys <- expand_surveys(set_den = c(1, 2) / 1000,
#'                           lengths_cap = c(100, 500),
#'                           ages_cap = c(5, 20))
#'
#' ## This call runs 500 simulations of 8 different surveys over the same
#' ## population, and then runs a stratified analysis and compares true vs
#' ## estimated values. It may take a while to run.
#' tests <- test_surveys(pop, surveys = surveys, keep_details = 1,
#'                       n_sims = 10, n_loops = 50, cores = 2)
#'
#' library(plotly)
#' tests$total_strat_error %>%
#'     filter(survey == 8, sim %in% 1:50) %>%
#'     group_by(sim) %>%
#'     plot_ly(x = ~year) %>%
#'     add_lines(y = ~I_hat, alpha = 0.5, name = "estimated") %>%
#'     add_lines(y = ~I, color = I("black"), name = "true") %>%
#'     layout(xaxis = list(title = "Year"),
#'            yaxis = list(title = "Abundance index"))
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
#'
#' }
#'
#' @export
#'
#' @import progress
#' @import doParallel
#' @import parallel
#' @import foreach
#'

test_surveys <- function(sim, surveys = expand_surveys(), keep_details = 1,
                         n_sims = 1, n_loops = 100, cores = 2, export_dir = NULL,
                         length_group = "inherit", alk_scale = "division", ...) {

  if (!is.null(export_dir)) {
    save(list = ls(all.names = TRUE),
         file = file.path(export_dir, "test_inputs.RData"))
  }
  .test_loop(sim = sim, surveys = surveys, n_sims = n_sims,
             n_loops = n_loops, cores = cores, export_dir = export_dir,
             complete = NULL, keep_details = keep_details,
             length_group = length_group, alk_scale = alk_scale, ...)

}


#' @export
#' @rdname test_surveys
resume_test <- function(export_dir = NULL) {
  load(file.path(export_dir, "test_inputs.RData"))
  load(file.path(export_dir, "complete.RData"))
  .test_loop(sim = sim, surveys = surveys, n_sims = n_sims,
             n_loops = n_loops, cores = cores, export_dir = export_dir,
             complete = complete, keep_details = keep_details,
             length_group = length_group, alk_scale = alk_scale, ...)
}



## Helper function for test_surveys and resume_test
.test_loop <- function(sim = NULL, surveys = NULL, n_sims = NULL,
                       n_loops = NULL, cores = NULL, export_dir = NULL,
                       complete = NULL, keep_details = NULL,
                       length_group = NULL, alk_scale = NULL, ...) {

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

    cat("\nRunning simulations...\n")
    pb <- progress::progress_bar$new(
      format = "Survey :current of :total [:bar] :percent in :elapsed (eta: :eta)",
      total = nrow(surveys), clear = FALSE, show_after = 0, width = 100)
    invisible(pb$tick(min(incomplete) - 1))

    for (i in incomplete) {

      cl <- makeCluster(cores) # use parallel computation
      registerDoParallel(cl)
      loop_error <- foreach(j = seq(n_loops), .packages = "SimSurvey") %dopar% {
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
      stopCluster(cl) # stop parallel process

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

      pb$tick()

    }

  }

  ## Compile results
  cat("\nCompiling results...\n")
  if (!is.null(export_dir)) {
    pb <- progress::progress_bar$new(
      format = "Survey :current of :total [:bar] :percent in :elapsed (eta: :eta)",
      total = nrow(surveys), clear = FALSE, show_after = 0, width = 100)
    for (i in seq(nrow(surveys))) {
      samp_totals[[i]] <- data.table::fread(file.path(export_dir, "samp_totals", paste0(i, ".csv")),
                                            verbose = FALSE, showProgress = FALSE)
      total_strat_error[[i]] <- data.table::fread(file.path(export_dir, "total_strat_error", paste0(i, ".csv")),
                                                  verbose = FALSE, showProgress = FALSE)
      length_strat_error[[i]] <- data.table::fread(file.path(export_dir, "length_strat_error", paste0(i, ".csv")),
                                                   verbose = FALSE, showProgress = FALSE)
      age_strat_error[[i]] <- data.table::fread(file.path(export_dir, "age_strat_error", paste0(i, ".csv")),
                                                verbose = FALSE, showProgress = FALSE)
      pb$tick()
    }
  }
  sim$surveys <- data.table::data.table(surveys)
  sim$samp_totals <- data.table::rbindlist(samp_totals); rm(samp_totals)
  sim$total_strat_error <- data.table::rbindlist(total_strat_error); rm(total_strat_error)
  sim$length_strat_error <- data.table::rbindlist(length_strat_error); rm(length_strat_error)
  sim$age_strat_error <- data.table::rbindlist(age_strat_error); rm(age_strat_error)

  ## Calculate some stats (limit to RMSE to limit size)
  sim$total_strat_error_stats <- sim$total_strat_error[, list(RMSE = sqrt(mean(error ^ 2))), by = "survey"]
  sim$length_strat_error_stats <- sim$length_strat_error[, list(RMSE = sqrt(mean(error ^ 2))), by = "survey"]
  sim$age_strat_error_stats <- sim$age_strat_error[, list(RMSE = sqrt(mean(error ^ 2))), by = "survey"]

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




