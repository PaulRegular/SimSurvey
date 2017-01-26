## Required libraries for data processing
library(sp)
library(rgdal)
library(dismo)
library(raster)
library(rgeos)
library(ggplot2)

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
strat_polys <- readOGR("data-raw/DFO_NL_survey_strat/2HJ3KLNOP_Strata_Polygons_WGS84.shp",
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
can <- raster::getData("GADM", country = "CAN", level = 1, path = "data-raw/GADM_data")
nl <- can[can$NAME_1 == "Newfoundland and Labrador", ]
spm <- raster::getData("GADM", country = "SPM", level = 1, path = "data-raw/GADM_data")
row.names(nl) <- paste("NL", row.names(nl), sep = "_") # provide unique ID's before attempting rbind
row.names(spm) <- paste("SPM", row.names(spm), sep = "_")
nl <- rbind(nl, spm)
nl_utm <- spTransform(nl, CRS(utm_proj))
nl_utm <- raster::crop(nl_utm, survey_extent)
nl_utm <- rgeos::gSimplify(nl_utm, tol = 0.2) # simplify the land object
plot(nl_utm)

## Import GEBCO derived bathymetry data
## "NAFO_GEBCO_data.nc" file was downloaded from http://www.gebco.net/.
##  Only provide derived data in data-raw folder for licensing reasons
# bathy <- raster("data-raw/GEBCO_derived_bathy/NAFO_GEBCO_data.nc")
# bathy <- raster::crop(bathy, extent(-65, -45, 42, 50)) # coarse crop
# bathy_utm <- raster::projectRaster(bathy, crs = utm_proj)
# r <- raster(survey_extent, res = c(1, 1), crs = utm_proj)
# bathy_utm <- resample(bathy_utm, r, method = "bilinear") # interpolate to 1 x 1 km utm grid
# save(bathy_utm, file = "data-raw/GEBCO_derived_bathy/GEBCO_derived_bathy.Rdata")
load("data-raw/GEBCO_derived_bathy/GEBCO_derived_bathy.Rdata")

## Import 3PS survey units
survey_units <- read.table("data-raw/DFO_NL_survey_units/POS3P.DAT", header = FALSE)
names(survey_units) <- c("division", "strat", "unit_num", "lat", "lon")
survey_units <- survey_units[survey_units$strat %in% index_strata, ]
survey_units$lat <- llconvert(survey_units$lat)
survey_units$lon <- -llconvert(survey_units$lon)
coordinates(survey_units) <- ~ lon + lat
proj4string(survey_units) <- proj4string(strat_polys)
survey_units_utm <- spTransform(survey_units, crs(utm_proj))

## Clean-up some unused objects
rm(strat_polys, nl, survey_units)


## Identify the survey region and add a buffer zone
## The buffer zone will minimize or get rid of edge effects in the survey
## area of the simulation + improve distance through water calculations
## (buffer zone may not be necessary for some cases)
survey_region <- raster::buffer(strat_polys_utm, width = 0.5) # km
buffer_zone <- raster::buffer(survey_region, width = 25)
<<<<<<< HEAD
=======
  samp_area <- raster::area(survey_region)

## Generate a regular grid over the survey & buffer zones
samp_den <- length(survey_units_utm) / samp_area
buffer_area <- raster::area(buffer_zone)
nunits <- round(buffer_area * samp_den)
survey_grid <- spsample(buffer_zone, n = nunits, type = "regular")
survey_grid <- as(as(survey_grid, "SpatialPixels"), "SpatialPolygons")
>>>>>>> 7d1e093a7031205072dce74bb520cfaf3cbb4b81

## Erase the land from the buffer_zone and clip the survey_grid
buffer_zone <- raster::erase(buffer_zone, nl_utm)
plot(buffer_zone, col = "grey") # clean up isolated polygon
# locator() # used this function to generate values below
p <- cbind(c(731.4171, 740.1683, 762.3381, 750.0864, 731.4171),
           c(5289.713, 5296.714, 5270.461, 5262.876, 5289.713))
p <- SpatialPolygons(list(Polygons(list(Polygon(p)), "1")),
                     pO = 1L, proj4string = CRS(proj4string(buffer_zone)))
# p <- raster::drawPoly() # easier (less repeatable) option
buffer_zone <- raster::erase(buffer_zone, p)
plot(buffer_zone, col = "grey")
<<<<<<< HEAD

## Generate a regular grid over the area, then clip to buffer_zone
samp_area <- raster::area(survey_region)
samp_den <- length(survey_units_utm) / samp_area
extent_poly <- as(survey_extent, "SpatialPolygons")
proj4string(extent_poly) <- proj4string(survey_region)
extent_area <- raster::area(extent_poly)
nunits <- round(extent_area * samp_den)
survey_grid <- spsample(extent_poly, n = nunits, type = "regular")
survey_grid <- as(as(survey_grid, "SpatialPixels"), "SpatialPolygons")
=======
  >>>>>>> 7d1e093a7031205072dce74bb520cfaf3cbb4b81
survey_grid <- raster::intersect(survey_grid, buffer_zone)
plot(survey_grid)

## Extract the strat each centroid lies over
coords <- data.frame(coordinates(survey_grid))
coordinates(coords) <- coords
proj4string(coords) <- proj4string(survey_grid)
data <- over(coords, strat_polys_utm)
data <- data[, c("DIV", "STRAT")]
names(data) <- c("division", "strat")
survey_grid <- SpatialPolygonsDataFrame(survey_grid, data = data)
plot(survey_grid[!is.na(survey_grid$strat), ])
cols <- rainbow(length(unique(survey_grid$strat)))
plot(survey_grid, col = cols[factor(survey_grid$strat)]) # coarse

## Extract coords and area, and calculate mean depth
xy <- coordinates(survey_grid) # these centroids are essentially the actual survey units
survey_grid$easting <- xy[, 1] # + the buffer centroids
survey_grid$northing <- xy[, 2]
survey_grid$area <- raster::area(survey_grid)
survey_grid$depth <- raster::extract(bathy_utm, survey_grid, fun = mean)
survey_grid$id <- row.names(survey_grid)

## Quick plot of map data for a visual check of the grid
png("data-raw/survey_grid.png", height = 20, width = 20, units = "in", res = 600)
plot(bathy_utm, col = colorRampPalette(c("steelblue", "white"))(100))
plot(nl_utm, add = TRUE, col = "grey", border = "grey50")
plot(strat_polys_utm, add = TRUE, border = "steelblue")
plot(survey_units_utm, add = TRUE, pch = 16, cex = 0.1, col = "red")
plot(survey_grid, add = TRUE, border = "steelblue", lwd = 0.1)
box()
dev.off()

## Generate another quick plot showing depth associated with each cell
# spplot(survey_grid[, "depth"]) # SLOW, use ggplot2
survey_grid_dat <- fortify(survey_grid)
survey_grid_dat <- merge(survey_grid_dat, survey_grid@data, by = "id")
ggplot(survey_grid_dat) +
  geom_polygon(aes(x = long, y = lat, group = group, fill = depth), colour = "steelblue")

## Simplify identifier by providing a unique string
survey_grid$cell <- seq_along(survey_grid)
row.names(survey_grid) <- as.character(survey_grid$cell)

## Finalize Rdata object - order and rename cols, and round double values
survey_grid@data <- survey_grid@data[, c("cell", "division", "strat",
                                         "easting", "northing", "area", "depth")]
survey_grid@data[] <- lapply(names(survey_grid@data),
                             function(nm) {
                               if(nm %in% c("easting", "northing", "area", "depth")) {
                                 round(survey_grid@data[[nm]], 2)
                               } else {
                                 survey_grid@data[[nm]]
                               }
                             }) # two digits should suffice for the double columns
save(survey_grid, file = "data/survey_grid.rda")


## Also export land and bathy object (optional, for plotting)
land <- nl_utm
save(land, file = "data/land.rda")
bathy <- bathy_utm
save(bathy, file = "data/bathy.rda")


