## Required libraries for data processing
library("sp")
library("rgdal")
library("raster")

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

## Build survey grid - approximate the resolution of the actual survey
plot(strat_polys_utm)
plot(survey_units_utm, add = TRUE, pch = ".", col = "red")
nunits <- length(survey_units_utm)
survey_grid <- spsample(strat_polys_utm, n = nunits, type = "regular")
plot(survey_grid, add = TRUE, col = "blue", pch = ".")
survey_grid <- as(survey_grid, "SpatialPixels")
dat <- data.frame(cell = seq.int(length(survey_grid)))
survey_grid <- SpatialPixelsDataFrame(survey_grid, dat)
survey_grid <- raster(survey_grid)
names(survey_grid) <- "cell"

## Add division, strat and depth
survey_grid$division <- rasterize(strat_polys_utm, survey_grid, "DIV")
survey_grid$strat <- rasterize(strat_polys_utm, survey_grid, "STRAT")
survey_grid$depth <- resample(bathy_utm, survey_grid, method = "bilinear")
values(survey_grid$depth)[is.na(values(survey_grid$cell))] <- NA
plot(survey_grid)

## Make sure the number of cells are equal across the stack
lapply(names(survey_grid), function(nm) sum(!is.na(values(survey_grid[[nm]]))))

## Export survey_grid + land bathy object (optional, for plotting)
save(survey_grid, file = "data/survey_grid.rda")
land <- nl_utm
save(land, file = "data/land.rda")
bathy <- bathy_utm
save(bathy, file = "data/bathy.rda")


