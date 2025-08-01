#' Make a flexdashboard for visualizing the simulation
#'
#' Launches an interactive flexdashboard to visualize simulation outputs.
#' Assumes the working directory is the root project directory.
#'
#' @param sim An object produced by [`sim_abundance()`], [`sim_distribution()`],
#' [`sim_survey()`], or [`test_surveys()`].
#' @param ... Additional arguments passed to [`rmarkdown::run()`].
#'
#' @return No return value. This function launches an interactive dashboard in the Viewer pane or browser.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   pop <- sim_abundance(ages = 1:20, years = 1:20)
#'   vis_sim(pop)
#'
#'   dist <- sim_distribution(pop, grid = make_grid(res = c(10, 10)))
#'   vis_sim(dist)
#'
#'   # Single survey
#'   survey <- sim_survey(dist, n_sims = 5)
#'   vis_sim(survey)
#'
#'   # Multiple survey designs
#'   surveys <- expand_surveys(set_den = c(1, 2) / 1000,
#'                             lengths_cap = c(100, 500),
#'                             ages_cap = c(5, 20))
#'
#'   tests <- test_surveys(dist, surveys = surveys, keep_details = 1,
#'                         n_sims = 5, n_loops = 5, cores = 1)
#'   vis_sim(tests)
#' }
#' }
#'
#' @export

vis_sim <- function(sim, ...) {

  pkg <- c("rmarkdown", "shiny", "flexdashboard", "crosstalk", "plotly", "viridis")
  for (p in pkg) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(paste(p, "is needed for vis_fit to work. Please install it."), call. = FALSE)
    }
  }

  rmd_env <- new.env()
  rmd_env$sim <- sim
  rmarkdown::run(file = system.file("rmd", "vis_sim.Rmd", package = "SimSurvey"),
                 render_args = list(envir = rmd_env), ...)

}

