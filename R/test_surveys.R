

#' Test sampling design of multiple surveys
#'
#' @description   This function allows a series of sampling design settings to
#' be set and tested on the simulated population. All combinations of the
#' supplied settings (\code{set_den}, \code{lengths_cap}, \code{ages_cap}) are tested.
#'
#' @param sim               Simulation from \code{\link{sim_distribution}}.
#' @param n_sims            Number of times to simulate a survey over the simulated population.
#'                          Requesting a large number of simulations here may max out your RAM.
#' @param n_loops           Number of times to run the \code{\link{sim_survey}} function. Total
#'                          simulations run will be the product of \code{n_sims} and \code{n_loops}
#'                          arguments. Low numbers of \code{n_sims} and high numbers of \code{n_loops}
#'                          will be easier on RAM, but slower.
#' @param set_den           Vector of set densities (number of sets per [grid unit] squared)
#' @param lengths_cap       Vector of maximum number of lengths measured per set
#' @param ages_cap          Vector of maximum number of otoliths to collect per length group
#'                          per division per year
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
#'

test_surveys <- function(sim, n_sims = 10, n_loops = 10,
                         set_den = c(0.3, 0.5, 0.8, 1, 2, 3, 6, 9) / 1000,
                         lengths_cap = c(2, 3, 5, 8, 10, 20, 30, 60, 90, 100, 200, 400, 600, 1000),
                         ages_cap = c(2, 3, 5, 8, 10, 20, 30, 60),
                         ...) {

  surveys <- expand.grid(set_den = set_den, lengths_cap = lengths_cap, ages_cap = ages_cap)
  surveys$survey <- seq(nrow(surveys))

  survey_error <- vector("list", nrow(surveys))

  pb <- progress::progress_bar$new(
    format = "Running [:bar] :percent in :elapsed (eta: :eta)",
    total = nrow(surveys), clear = FALSE, width = 100,
    show_after = 0)
  invisible(pb$tick(0))

  for (i in surveys$survey) {

    loop_error <- vector("list", n_loops)
    sim_counter <- 0

    for (j in seq(n_loops)) {
      res <- sim_survey(sim, n_sims = n_sims,
                        set_den = surveys$set_den[i],
                        lengths_cap = surveys$lengths_cap[i],
                        ages_cap = surveys$ages_cap[i],
                        ...) %>%
        run_strat() %>% strat_error()
      total_strat_error <- res$total_strat_error
      total_strat_error$sim <- total_strat_error$sim + sim_counter
      age_strat_error <- res$age_strat_error
      age_strat_error$sim <- age_strat_error$sim + sim_counter
      sim_counter <- sim_counter + n_sims
      loop_error[[j]]$total_strat_error <- total_strat_error
      loop_error[[j]]$age_strat_error <- age_strat_error
    }

    total_strat_error <- data.table::rbindlist(lapply(loop_error, `[[`, "total_strat_error"))
    age_strat_error <- data.table::rbindlist(lapply(loop_error, `[[`, "age_strat_error"))
    total_strat_error$survey <- i
    age_strat_error$survey <- i
    survey_error[[i]]$total_strat_error <- total_strat_error
    survey_error[[i]]$age_strat_error <- age_strat_error
    survey_error[[i]]$total_strat_error_stats <- c(survey = i, error_stats(total_strat_error$error))
    survey_error[[i]]$age_strat_error_stats <- c(survey = i, error_stats(age_strat_error$error))

    pb$tick()

  }

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

