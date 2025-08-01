
#' Helper function to generate precision matrix Q for simulation
#'
#' `r lifecycle::badge("experimental")`
#'
#' Creates a precision matrix from a mesh created by **R-INLA**, given a specified
#' decorrelation range. This function supports both the standard SPDE approach and
#' the barrier model version. It serves a similar purpose to `.sp_covar()`.
#'
#' @param mesh The mesh created by R-INLA representing the spatial area.
#' @param barrier.triangles A list of mesh triangles that represent the barrier (used only for the barrier model).
#' @param range The decorrelation range.
#' @param range_fraction Fraction used to adjust the range across the "land" (barrier) regions. Only used in the barrier model.
#' @param sigma_u The overall variance of the spatial process.
#' @param model Either `"spde"` for the standard SPDE model or `"barrier"` for the barrier model.
#'
#' @return A sparse precision matrix `Q` of class `dgTMatrix`.
#'
#' @keywords internal
.Q <- function(mesh,barrier.triangles, range = 50, range_fraction = 0.2,
               sigma_u = 1, model = "spde"){

  Q <-  switch(model,
               spde = {
                 ##Priors aren't used...
                 ##I just used this so precision uses range...
                 i = INLA::inla.spde2.pcmatern(mesh, prior.range = c(1, 0.1), prior.sigma = c(1, 0.1))
                 INLA::inla.spde2.precision(i, theta = c(log(range), log(sigma_u)))
               },
               barrier = {

                 # ##Switch out from using inla.barrier.fem and :::, not really needed here anyways
                 # barrier.model <- INLA::inla.barrier.pcmatern(mesh,barrier.triangles=barrier.triangles)
                 # ##Yes, theta really is supposed to be the reverse of the one above...
                 # INLA::inla.rgeneric.q(barrier.model,"Q",theta=c(log(sigma_u),log(range)))

                 ## Using 'new' approach described in https://eliaskrainski.github.io/INLAspacetime/articles/web/barrierExample.html
                 barrier.model <- INLAspacetime::mesh2fem.barrier(mesh, barrier.triangles)
                 INLA::inla.barrier.q(barrier.model, ranges = c(range, range * range_fraction), sigma = sigma_u)

               },
               stop("wrong or no specification of covariance model"))
  Q
}

#' Simulate age-year-space covariance using SPDE approach
#'
#' `r lifecycle::badge("experimental")`
#'
#' Returns a function for use inside [`sim_distribution()`] to generate the error term.
#'
#' @param sd Variance of the process (can be age-specific).
#' @param range Decorrelation range.
#' @param model Either `"barrier"` or `"spde"`; determines how the precision matrix `Q` is generated.
#' @param phi_age Autocorrelation through ages. Can be a single value or a vector the same length as `ages`.
#' @param phi_year Autocorrelation through years. Can be a single value or a vector the same length as `years`.
#' @param group_ages Ages to group together for shared space-age-year variance.
#' @param group_years Years to group together for shared space-age-year variance.
#' @param mesh The mesh used to generate the precision matrix.
#' @param barrier.triangles The set of mesh triangles that define the barrier (used only in the barrier model).
#'
#' @return A function that can be passed to [`sim_distribution()`] as the `ays_covar` argument.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA")) {
#'
#'   # Make a grid
#'   my_grid <- make_grid(res = c(10, 10))
#'
#'   # Make a mesh based on the grid
#'   my_mesh <- make_mesh(my_grid)
#'
#'   # Simulate and plot
#'   sim <- sim_abundance(ages = 1:10, years = 1:10) |>
#'     sim_distribution(
#'       grid = my_grid,
#'       ays_covar = sim_ays_covar_spde(
#'         phi_age = 0.8,
#'         phi_year = 0.1,
#'         model = "spde",
#'         mesh = my_mesh
#'       ),
#'       depth_par = sim_parabola(mu = 200, sigma = 50)
#'     )
#'
#'   plot_distribution(sim, ages = 1:5, years = 1:5, type = "heatmap")
#' }
#' }
#'
#' @export

sim_ays_covar_spde <- function(sd = 2.8,
                               range = 300,
                               model = "spde",
                               phi_age = 0.5,
                               phi_year = 0.9,
                               group_ages = 5:20,
                               group_years = NULL,
                               mesh,
                               barrier.triangles) {

  function(x = NULL, ages = NULL, years = NULL, cells = NULL){

    for (pkg in c("INLA")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste(pkg, "is needed for make_mesh to work. Please install it."), call. = FALSE)
      }
    }

    na <- length(ages)
    ny <- length(years)
    nc <- length(cells)
    if (length(sd) == 1) {
      sd <- rep(sd, na)
    } else {
      if (length(sd) != na) {
        stop("The number of sd values supplied != number of ages.")
      }
    }
    if (length(phi_age) == 1) {
      phi_age <- rep(phi_age, na)
    } else {
      if (length(phi_age) != na) {
        stop("The number of phi_age values supplied != number of ages.")
      }
    }
    if (length(phi_year) == 1) {
      phi_year <- rep(phi_year, ny)
    } else {
      if (length(phi_year) != ny) {
        stop("The number of phi_year values supplied != number of years.")
      }
    }
    age_map <- as.character(ages)
    if (!is.null(group_ages)) {
      age_map[ages %in% group_ages] <- paste(range(group_ages), collapse = ":")
    }
    year_map <- as.character(years)
    if (!is.null(group_years)) {
      year_map[years %in% group_years] <- paste(range(group_years), collapse = ":")
    }

    pc_age <- sqrt(1 - phi_age ^ 2)
    pc_year <- sqrt(1 - phi_year ^ 2)

    Q <- .Q(mesh,barrier.triangles,range,model=model)
    ##A matrix to relate the cell locations to the mesh locations
    A <- INLA::inla.spde.make.A(mesh,as.matrix(x))
    ##Draw samples from GF defined by Q and project it to the cell locs
    u.data <- as.matrix(A%*%INLA::inla.qsample(na*ny,Q=Q))
    ##Fill up the error array with the samples
    E <- array(t(u.data),dim=c(na,ny,nc),
               dimnames=list(age = ages,year = years,cell=cells))


    ## Transform them to the distribution defined in Append. 1 of
    ## the Sim Survey paper
    for(j in seq_along(years)){
      for(i in seq_along(ages)){
        if((i == 1) & (j == 1)){
          m <- 0
          s <- sd[i]/(pc_age[i]*pc_year[j])
          E[i,j,] = E[i,j,]*s + m
        }
        if((i > 1) & (j == 1)){
          if(age_map[i] == age_map[i-1]){
            E[i,j,] <- E[i-1,j,]
          } else{
            m <- phi_age[i] * E[i-1,j,]
            s <- sd[i]/pc_year[j]
            E[i,j,] = E[i,j,]*s + m
          }

        }
        if((i == 1) & (j > 1)){
          if(year_map[j] == year_map[j-1]){
            E[i,j,] <- E[i,j-1,]
          } else{
            m <- phi_year[j] * E[i,j-1,]
            s <- sd[i]/pc_age[i]
            E[i,j,] <- E[i,j,]*s +m

          }
        }
        if((i > 1) & (j > 1)){
          if (age_map[i] == age_map[i - 1]) {
            E[i, j, ] <- E[i - 1, j, ]
          }
          if (year_map[j] == year_map[j - 1]) {
            E[i, j, ] <- E[i, j - 1, ]
          }
          if ((age_map[i] != age_map[i - 1]) & (year_map[j] != year_map[j - 1])) {
            m <- phi_year[j] * E[i, j - 1, ] + phi_age[i] * (E[i - 1, j, ] - phi_year[j] * E[i - 1, j - 1, ])
            s <- sd[i]
            E[i,j,] <- E[i,j,]*s + m

          }
        }
      }
    }
    E
  }

}


#' Make an R-INLA mesh based on a grid
#'
#' This function creates a mesh based on a given grid. While mesh construction
#' and validation should ideally be done manually, this function provides a
#' convenient default interface between a grid and `inla.mesh.2d`. It is designed
#' to support usage with [`sim_ays_covar_spde()`], and the default parameters are
#' tuned for use with the default grid setup.
#'
#' @param grid A grid object to generate a mesh from.
#' @param max.edge The maximum allowed triangle edge length. One or two numeric values.
#' Passed to `inla.mesh.2d`.
#' @param bound.outer Optional outer extension value passed to `offset`.
#' @param offset Automatic extension distance used by `inla.mesh.2d`.
#' @param cutoff Minimum distance allowed between mesh points.
#' @param ... Additional options passed to `inla.mesh.2d`.
#'
#' @return An object of class `inla.mesh`.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA")) {
#'   basic_mesh <- make_mesh()
#'   plot(basic_mesh)
#' }
#' }
#'
#' @export

make_mesh <- function(grid = make_grid(),
                      max.edge = 50,
                      bound.outer = 150,
                      cutoff = 10,
                      offset = c(max.edge, bound.outer),
                      ...) {

  for (pkg in c("INLA")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste(pkg, "is needed for make_mesh to work. Please install it."), call. = FALSE)
    }
  }

  gridPoints <- data.frame(grid)
  locs <- as.matrix(gridPoints[,1:2])
  mesh <- INLA::inla.mesh.2d(locs, offset = offset, max.edge = max.edge, cutoff = cutoff, ...)
  mesh

}





