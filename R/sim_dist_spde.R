
#' Helper function to generate precision matrix Q for simulation
#'
#' This creates the precision matrix from a mesh created by R-INLA with a
#' specified range. This is currently set up to support the standard spde
#' approach of precision matrices and the barrier model version. Similar
#' in purpose to .sp_covar.
#'
#' @param mesh The mesh created by R-INLA representing spatial area
#' @param barrier.triangles the list of triangles of the mesh in the barrier, only used for barrier model
#' @param range decorrelation range
#' @param range_fraction the fraction that sets the range over the "land"
#' parts of the barrier model. Only used with the barrier model
#' @param sigma_u the overall variance of the spatial process
#' @param model 'spde' for the SPDE approach, 'barrier' for barrier approach
#'
#' @return Q a sparse precision matrix of type dgTMatrix
#' @keywords internal
#'
.Q <- function(mesh,barrier.triangles,range=50,range_fraction=0.2,
               sigma_u = 1,model="spde"){
  Q <-  switch(model,
               spde = {
                 ##Priors aren't used...
                 ##I just used this so precision uses range...
                 i = INLA::inla.spde2.pcmatern(mesh,prior.range=
                                                 c(1,0.1),
                                               prior.sigma=c(1,0.1))
                 INLA::inla.spde2.precision(i,theta=c(log(range),log(sigma_u)))
               },
               barrier = {
                 ##Switch out from using inla.barrier.fem and :::, not really needed here anyways
                 barrier.model <- INLA::inla.barrier.pcmatern(mesh,barrier.triangles=barrier.triangles)
                 ##Yes, theta really is supposed to be the reverse of the one above...
                 INLA::inla.rgeneric.q(barrier.model,"Q",theta=c(log(sigma_u),log(range)))
               },
               stop("wrong or no specification of covariance model"))
  Q
}

#' Simulate age-year-space covariance using SPDE approach
#'
#' Returns a function to use inside \code{\link{sim_distribution}} to
#' generate the error term.
#'
#' @param sd Variance (can be age specific)
#' @param range Decorrelation range
#' @param model String indicating "barrier" or "spde" to generate Q with
#' @param phi_age Defines autocorrelation through ages. Can be one value or
#'        a vector of the same length as ages.
#' @param phi_year Defines autocorrelation through years. Can be one value
#'        or a vector of the same length as years.
#' @param group_ages Make space-age-year variance equal across these ages
#' @param group_years Make space-age-year variance equal across these years
#' @param mesh The mesh used to generate the precision matrix
#' @param barrier.triangles the set of triangles in the barrier of the mesh
#' for the barrier model
#'
#' @return Returns a function for use in \code{\link{sim_distribution}}.
#'
#' @examples
#'
#' ##SPDE Approach
#'
#' \donttest{
#'
#' ## Make a grid
#' my_grid <- make_grid(res = c(10,10))
#'
#' ## Make a mesh based off it
#'
#' my_mesh <- make_mesh(my_grid)
#' sim <- sim_abundance(ages = 1:10, years = 1:10) %>%
#'          sim_distribution(grid = my_grid,
#'                           ays_covar = sim_ays_covar_spde(phi_age = 0.8,
#'                                                          phi_year = 0.1,
#'                                                          model = "spde",
#'                                                          mesh = my_mesh),
#'                           depth_par = sim_parabola(mu = 200,
#'                                                    sigma = 50))
#' plot_distribution(sim,ages = 1:5, years = 1:5, type = "heatmap")
#'
#' }
#'
#' @export
#'

sim_ays_covar_spde <- function(sd = 2.8,
                               range = 300,
                               model = "barrier",
                               phi_age = 0.5,
                               phi_year = 0.9,
                               group_ages = 5:20,
                               group_years = NULL,
                               mesh,
                               barrier.triangles) {

  function(x = NULL, ages = NULL, years = NULL, cells = NULL){

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


#' Make an R-INLA mesh based off a grid
#'
#' This will make a mesh based off a given grid. Ideally the mesh
#' construction and validation should be done by hand, but this exists
#' for convenience. Meshes are used for sim_ays_covar_spde. The defaults
#' are designed for the default grid. Just a basic interface between the grid and inla.mesh.2d.
#'
#' @param grid grid object to make a mesh of
#' @param max.edge The largest allowed triangle edge length. One or two values. This is passed to inla.mesh.2d
#' @param bound.outer The optional outer extension value given to offset.
#' @param offset The automatic extension distance given to inla.mesh.2d
#' @param cutoff Minimum distance allowed between points
#' @param ... Other options to pass to inla.mesh.2d
#'
#' @examples
#'
#' \donttest{
#'
#' basic_mesh <- make_mesh()
#' plot(basic_mesh)
#'
#' }
#'
#' @export
#'

make_mesh <- function(grid = make_grid(),
                      max.edge = 50,
                      bound.outer = 150,
                      cutoff = 10,
                      offset = c(max.edge, bound.outer),
                      ...) {

  for (pkg in c("rgdal", "INLA")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste(pkg, "is needed for make_mesh to work. Please install it."), call. = FALSE)
    }
  }

  gridPoints <- data.frame(grid)
  locs <- as.matrix(gridPoints[,1:2])
  mesh <- INLA::inla.mesh.2d(locs, offset = offset, max.edge = max.edge, cutoff = cutoff, ...)
  mesh

}





