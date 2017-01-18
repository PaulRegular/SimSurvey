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
index.strata <- c(293:300, 306:326, 705:708, 711:716, 779:783)

## UTM projection for Newfoundland
utm.proj <- "+proj=utm +zone=21 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

## Import 3Ps strata shapefile
strat.polys <- readOGR("data-raw/DFO_NL_survey_strat/2HJ3KLNOP_Strata_Polygons_WGS84.shp", layer = "2HJ3KLNOP_Strata_Polygons_WGS84")
strat.polys <- strat.polys[strat.polys$DIV == "3P", ] # subset to 3P
strat.polys <- strat.polys[strat.polys$STRAT %in% index.strata, ] # 3Ps index strata
strat.polys.utm <- spTransform(strat.polys, utm.proj)

## Import GEBCO derived bathymetry data
## "NAFO_GEBCO_data.nc" file was downloaded from http://www.gebco.net/.
##  Only provide derived data in data-raw folder for licensing reasons
# bathy <- raster("data-raw/GEBCO_derived_bathy/NAFO_GEBCO_data.nc")
# bathy <- raster::crop(bathy, extent(strat.polys))
# bathy.utm <- raster::projectRaster(bathy, crs = utm.proj)
# r <- raster(extent(bathy.utm), res = c(1, 1), crs = utm.proj)
# bathy.utm <- resample(bathy.utm, r, method = "bilinear") # interpolate to 1 x 1 km utm grid
# bathy <-  raster::projectRaster(bathy.utm, crs = proj4string(strat.polys))
# save(bathy, bathy.utm, file = "GEBCO_derived_bathy/GEBCO_derived_bathy.Rdata")
load("data-raw/GEBCO_derived_bathy/GEBCO_derived_bathy.Rdata")

## Import land data
can <- raster::getData("GADM", country = "CAN", level = 1, path = "data-raw/GADM_data")
nl <- can[can$NAME_1 == "Newfoundland and Labrador", ]
spm <- raster::getData("GADM", country = "SPM", level = 1, path = "data-raw/GADM_data")
row.names(nl) <- paste("NL", row.names(nl), sep = "_") # provide unique ID's before attempting rbind
row.names(spm) <- paste("SPM", row.names(spm), sep = "_")
nl <- rbind(nl, spm)
nl.utm <- spTransform(nl, CRS(utm.proj))

## Import 3PS survey units
survey.units <- read.table("data-raw/DFO_NL_survey_units/POS3P.DAT", header = FALSE)
names(survey.units) <- c("division", "strat", "unit.num", "lat", "lon")
survey.units <- survey.units[survey.units$strat %in% index.strata, ]
survey.units$lat <- llconvert(survey.units$lat)
survey.units$lon <- -llconvert(survey.units$lon)
coordinates(survey.units) <- ~ lon + lat
proj4string(survey.units) <- proj4string(strat.polys)
survey.units.utm <- spTransform(survey.units, crs(utm.proj))

## Clean-up one issue with the strat polygons in the loop below
plot(strat.polys.utm[strat.polys.utm$Strat_ID %in% c(295, 296), ],
     col = rgb(1, 0, 0, 0.5)) # overlap!
plot(survey.units.utm[survey.units.utm$strat %in% c(295, 296), ],
     add = TRUE, pch = 16, col = rgb(1, 0, 0, 0.5))

## Make an irregular grid using Voronoi tessellation
## Loop across each strata to impose proper bounds
simgrids <- vector("list", length(index.strata))
for(i in seq_along(index.strata)) {
  u <- survey.units.utm[survey.units.utm$strat == index.strata[i], ]
  s <- strat.polys.utm[strat.polys.utm$Strat_ID == index.strata[i], ]
  s <- raster::buffer(s, width = 0.00001) # add a small buffer to fix some issues
  if(index.strata[i] == 295) {
    s2 <- strat.polys.utm[strat.polys.utm$Strat_ID == 296, ]
    s2 <- raster::buffer(s2, width = 0.00001)
    s <- raster::erase(s, s2)
  }
  v <- dismo::voronoi(u, ext = extent(s))
  v <- raster::intersect(s, v)
  row.names(v) <- paste0(v$strat, "-", v$unit.num) # ensure unique row names are provided
  simgrids[[i]] <- v
  plot(v, main = index.strata[i])
  plot(s, add = TRUE, border = "red", lty = 3)
  plot(u, add = TRUE, pch = 16, cex = 0.5)
}
simgrid <- do.call(rbind, simgrids)
## There are some allignment issues here, but they shouldn't be an issue

## Quick plot of map data for a visual check of the grid
png("data-raw/simgrid.png", height = 20, width = 20, units = "in", res = 600)
plot(bathy.utm, col = colorRampPalette(c("steelblue", "white"))(100))
plot(nl.utm, add = TRUE, col = "grey", border = "grey50")
plot(strat.polys.utm, add = TRUE, border = "steelblue")
plot(survey.units.utm, add = TRUE, pch = 16, cex = 0.1, col = "red")
plot(simgrid, add = TRUE, border = "steelblue", lwd = 0.1)
box()
dev.off()

## Extract area and calculate mean depth
simgrid$area <- raster::area(simgrid)
simgrid$depth <- extract(bathy.utm, simgrid, fun = mean)
simgrid$id <- row.names(simgrid)

## Generate another quick plot showing depth associated with each cell
# spplot(simgrid[, "depth"]) # SLOW, use ggplot2
simgrid.dat <- fortify(simgrid)
simgrid.dat <- merge(simgrid.dat, simgrid@data, by = "id")
ggplot(simgrid.dat) +
  geom_polygon(aes(x = long, y = lat, group = group, fill = depth), colour = "steelblue")

## Finalize Rdata object by adding survey unit coordinates, then order and rename cols
simgrid@data <- merge(simgrid@data,  as.data.frame(survey.units.utm),
                          by = c("division", "strat", "unit.num"), all.x = TRUE)
simgrid@data <- simgrid@data[, c("id", "division", "strat", "unit.num",
                                         "lon", "lat", "area", "depth")]
names(simgrid@data) <- c("id", "division", "strat", "unit.num",
                             "easting", "northing", "area", "depth")
simgrid@data[] <- lapply(names(simgrid@data),
                             function(nm) {
                               if(nm %in% c("easting", "northing", "area", "depth")) {
                                 round(simgrid@data[[nm]], 2)
                               } else {
                                 simgrid@data[[nm]]
                               }
                             }) # two digits should suffice for the double columns
save(simgrid, file = "data/simgrid.Rdata")

