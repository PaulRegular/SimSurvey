#' Make a flexdashboard for visualizing population simulation
#'
#' Assumes the working directory is the project directory
#'
#' @param sim            Object produced by \code{\link{sim_abundance}} or \code{\link{sim_distribution}}
#' @param ...            Additional arguments to send to \link[rmarkdown]{run}
#'
#' @export
#'

vis_sim <- function(sim, ...) {

  pkg <- c("rmarkdown", "shiny", "flexdashboard", "plotly", "viridis")
  for (p in pkg) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(paste(p, "is needed for vis_sim to work. Please install it."), call. = FALSE)
    }
  }

  ## add a temp object to the global environment, otherwise rmarkdown
  ## can't see these data
  assign("temp", sim, envir = globalenv())
  rmarkdown::run(file = "inst/rmd/vis_sim.Rmd", ...)

}
