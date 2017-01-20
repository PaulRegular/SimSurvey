#' Sample survey simulation grid.
#'
#' A exemplar for the structure of a survey grid object to supply to the functions
#' in this package.
#'
#' @format A SpatialPolygonsDataFrame with 8 variables:
#' \describe{
#'   \item{cell}{Survey cell identifier}
#'   \item{division}{NAFO division}
#'   \item{strat}{Survey strata number}
#'   \item{easting}{Easting coordinates from UTM Zone 21, datum WGS84, units = km}
#'   \item{northing}{Northing coordinates from UTM Zone 21, datum WGS84, units = km}
#'   \item{area}{Area of Voronoi cell (polygon) the sampling unit lies within}
#'   \item{depth}{Mean depth of the waters under each Voronoi cell (polygon), units = m}
#' }
#'
#' For further details on how this file was created, see the data-raw folder for
#' this package
#'
"survey_grid"
