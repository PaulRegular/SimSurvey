#' Sample survey simulation grid.
#'
#' A exemplar for the structure of a survey grid object to supply to the functions
#' in this package.
#'
#' @format A stars object with 4 attributes:
#' \describe{
#'   \item{cell}{Survey cell identifier}
#'   \item{division}{NAFO division}
#'   \item{strat}{Survey strata number}
#'   \item{depth}{Mean depth of the waters under each cell, units = m}
#' }
#'
#' For further details on how this file was created, see the data-raw folder for
#' this package.
#'
"survey_grid"


#' Southern Newfoundland coastline
#'
#' @format A sf object (MULTIPOLYGON)
#'
#' Derived from global administrative boundaries data (http://gadm.org/) downloaded
#' using the \code{\link{getData}} function. Details provided in the
#' data-raw folder for this package.
#'
"land"


#' Southern Newfoundland bathymetry
#'
#' @format A stars object
#'
#' Derived from data downloaded from http://www.gebco.net/. Details provided in
#' the data-raw folder for this package.
#'
"bathy"

#' Sample survey meshes and related items
#'
#'  @format A list containing the R-INLA survey mesh, the set of triangles in the barrier and the barrier polygons for plotting
#'
#' An example of a mesh containing barrier information for use with
#' sim_ays_covar_spde. Also derived from global administrative boundaries
#' data (http://gadm.org). Details on creation provided in the data-raw
#' folder of this package in the survey_mesh.R file. Includes the set
#' of barrier triangles needed to use the barrier approach, barrier
#' polygons for plotting and the set of triangles in the barrier.
"survey_mesh"

#' Lite sample survey mesh and related items
#'
#' @format A list containing the same items as survey_mesh, but with fewer nodes to save on computational time
#'
"survey_lite_mesh"


