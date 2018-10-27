#' Make a flexdashboard for visualizing the simulation
#'
#' Assumes the working directory is the project directory
#'
#' @param sim            Object produced by \code{\link{sim_abundance}}, \code{\link{sim_distribution}},
#'                       \code{\link{sim_survey}} or \code{\link{test_surveys}}.
#' @param ...            Additional arguments to send to \link[rmarkdown]{run}
#'
#' @examples
#'
#' \dontrun{
#'
#' pop <- sim_abundance(ages = 1:20, years = 1:20)
#' vis_sim(pop)
#'
#' dist <- sim_distribution(pop, grid = make_grid(res = c(10, 10)))
#' vis_sim(dist)
#'
#' ## Run one survey design
#' survey <- sim_survey(dist, n_sims = 5)
#' vis_sim(survey)
#'
#' ## Run several survey designs and assess stratified estimates
#' surveys <- expand_surveys(set_den = c(1, 2) / 1000,
#'                           lengths_cap = c(100, 500),
#'                           ages_cap = c(5, 20))
#' tests <- test_surveys(dist, surveys = surveys, keep_details = 1,
#'                       n_sims = 10, n_loops = 50, cores = 2)
#' vis_sim(tests)
#'
#' }
#'
#' @export
#'

vis_sim <- function(sim, ...) {

  pkg <- c("rmarkdown", "shiny", "flexdashboard", "crosstalk", "plotly", "viridis")
  for (p in pkg) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(paste(p, "is needed for vis_fit to work. Please install it."), call. = FALSE)
    }
  }

  ## make a tmp object in the global environment for use by rmarkdown
  ## (rmarkdown::run likes to operate on objects in the global environment)
  tmp <- as.list(environment())
  assign("tmp", tmp, globalenv())

  ## keep most objects created in the rmd file local (i.e. not in the global)
  local({
    rmarkdown::run(file = "inst/rmd/vis_sim.Rmd")
  })

}

