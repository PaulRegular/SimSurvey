##Setup a mesh from survey_grid for use testing the barrier simulation
library("sp")
library("rgdal")
library("raster")
library("sf")

## Function for converting long and lat values in DAT files
llconvert <- function(x) {
  ifelse(as.numeric(substring(x, 1, 2)) != 0,
         as.numeric(substring(x,1,2)) +
           (as.numeric(substring(x,3,4)) +
              as.numeric(substring(x,5,5))/10)/60,
         NA)
}


## Index strata for 3Ps
index_strata <- c(293:300, 306:326, 705:708, 711:716, 779:783)

## UTM projection for Newfoundland
utm_proj <- "+proj=utm +zone=21 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

## Import 3Ps strata shapefile
strat_polys <- readOGR("./DFO_NL_survey_strat/2HJ3KLNOP_Strata_Polygons_WGS84.shp",
                       layer = "2HJ3KLNOP_Strata_Polygons_WGS84")
strat_polys <- strat_polys[strat_polys$DIV == "3P", ] # subset to 3P
strat_polys <- strat_polys[strat_polys$STRAT %in% index_strata, ] # 3Ps index strata
strat_polys_utm <- spTransform(strat_polys, utm_proj)

## Survey area extent + a buffer
survey_extent <- bbox(strat_polys_utm)
b <- 120 # km buffer
survey_extent[, "min"] <- survey_extent[, "min"] - b
survey_extent[, "max"] <- survey_extent[, "max"] + b
survey_extent <- extent(survey_extent)

## Import land data
can <- raster::getData("GADM", country = "CAN", level = 1, path = "./GADM_data")
nl <- can[can$NAME_1 == "Newfoundland and Labrador", ]
spm <- raster::getData("GADM", country = "SPM", level = 1, path = "./GADM_data")
row.names(nl) <- paste("NL", row.names(nl), sep = "_") # provide unique ID's before attempting rbind
row.names(spm) <- paste("SPM", row.names(spm), sep = "_")
nl <- rbind(nl, spm)
nl_utm <- spTransform(nl, CRS(utm_proj))
nl_utm <- raster::crop(nl_utm, survey_extent)
nl_utm <- rgeos::gSimplify(nl_utm, tol = 0.2) # simplify the land object
plot(nl_utm)


##Convert the strat_polys_utm to sf because that's what I know better!
strat_polys_utm <- st_as_sf(strat_polys_utm)

##Generate the bounding box of the 3Ps area
bbox3Ps <- st_bbox(strat_polys_utm)

SPMNew <- st_crop(st_as_sf(nl_utm),strat_polys_utm)
SPMNew <- st_union(SPMNew)
bb <- st_as_sfc(bbox3Ps)
diff <- st_difference(bb,SPMNew)

##Create the survey mesh, I wanted a similar number of nodes to the number
##of cells that appear in survey_grid. This might be too computationally
##heavy if used in a model.
library(INLA)
max.edge = 8
bound.outer = 150
survey_mesh <- inla.mesh.2d(boundary=as_Spatial(diff),
                     max.edge = c(1,5)*max.edge,
                     cutoff = 2,
                     offset = c(max.edge,bound.outer))

##Get the barrier triangles
water.tri <- inla.over_sp_mesh(as_Spatial(diff),y=survey_mesh,
                               type = "centroid",ignore.CRS = TRUE)
num.tri = length(survey_mesh$graph$tv[,1])
survey_barrier_tri = setdiff(1:num.tri,water.tri)
survey_barrier_poly = inla.barrier.polygon(survey_mesh,
                                           barrier.triangles = survey_barrier_tri)

##A lighter version of the mesh for comparison/faster computational time
##~1500 vs. 5800 nodes
max.edge = 30
bound.outer = 150
survey_mesh_lite <- inla.mesh.2d(boundary=as_Spatial(diff),
                                 max.edge = c(1,5)*max.edge,
                                 cutoff = 2,
                                 offset = c(max.edge,bound.outer))

##Get the barrier triangles
water.tri <- inla.over_sp_mesh(as_Spatial(diff),y=survey_mesh_lite,
                               type = "centroid",ignore.CRS = TRUE)
num.tri = length(survey_mesh_lite$graph$tv[,1])
survey_lite_tri = setdiff(1:num.tri,water.tri)
survey_lite_poly = inla.barrier.polygon(survey_mesh,
                                           barrier.triangles = survey_lite_tri)

survey_lite_mesh <- list(mesh=survey_mesh_lite,barrier_tri=survey_lite_tri,barrier_poly=survey_lite_poly)
survey_mesh <- list(mesh=survey_mesh,barrier_tri=survey_barrier_tri,barrier_poly=survey_barrier_poly)

##Save the necessary stuff for actually using this mesh and what not
save(survey_mesh,survey_lite_mesh,
     file = "../data/survey_mesh.rda",compress = "xz")
