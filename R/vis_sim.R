#' Make a flexdashboard for visualizing the simulation
#'
#' Assumes the working directory is the project directory
#'
#' @param sim            Object produced by \code{\link{sim_abundance}}, \code{\link{sim_distribution}},
#'                       \code{\link{sim_survey}} or \code{\link{test_surveys}}.
#' @param ...            Additional arguments to send to \link[rmarkdown]{run}
#'
#' @export
#'

vis_sim <- function(sim, ...) {

  pkg <- c("rmarkdown", "shiny", "flexdashboard", "plotly", "viridis")
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

