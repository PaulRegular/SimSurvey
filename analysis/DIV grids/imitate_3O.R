
## Explore sampling data in 3O from the RV survey and try to imitate
## the data using SimSurvey

## Survey grid -----------------------------------------------------------------

library(sf)
library(raster)
library(SimSurvey)

## UTM projection for Newfoundland
utm_proj <- "+proj=utm +zone=21 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

## Import 3LNO strata shapefile
strat_polys <- st_read("data-raw/DFO_NL_survey_strat/2HJ3KLNOP_Strata_Polygons_WGS84.shp",
                       layer = "2HJ3KLNOP_Strata_Polygons_WGS84")
strat_polys <- strat_polys[strat_polys$DIV %in% c("3O"), ] # subset to 3O
strat_polys_utm <- st_transform(strat_polys, utm_proj)

## Import bathy data
strat_bathy <- raster("data-raw/GEBCO_derived_bathy/2HJ3KLNOP_GEBCO_bathy.nc")
strat_bathy <- raster::mask(strat_bathy, strat_polys) # cells within strata
# strat_bathy_utm <- raster::projectRaster(strat_bathy, crs = utm_proj)

## Depth range in the survey area
range(strat_bathy[], na.rm = TRUE)
# -2102    -48

## Number of strata
l_strata <- length(unique(strat_polys$STRAT))
l_strata
# 36

## Area of the strata NEEDS TO BE FIXED
range(sf::st_area(strat_polys_utm))
# Units: [km^2]
# 33.82768 8871.60986

## Survey area
survey_area <- sum(sf::st_area(strat_polys_utm))
survey_area
# 67958.99 [km^2]

## For make_grid arguments

#x_y_range
sqrt(survey_area)/2
# 130

## Mean area
# mean(sf::st_area(strat_polys_utm)) ## doesn't work because it includes STRATA with same names
survey_area/l_strata
# 1904.106 [km^2]

## TODO: Modify make_grid arguments to create a similar survey area



grid <- make_grid(x_range = c(-130, 130),
                  y_range = c(-130, 130),
                  res = c(3.5, 3.5),
                  shelf_depth = 200,
                  shelf_width = 100,
                  depth_range = c(0, 2100),
                  n_div = 1,
                  strat_breaks = seq(0, 1000, by = 20),
                  strat_splits = 1)
prod(res(grid)) * ncell(grid)
# 67081
length(unique(grid$strat))

mean(table(values(grid$strat)) * res(grid))
#length must be even number or will come up with error

range(values(grid$depth), na.rm = TRUE)

plot(grid)
plot(grid$depth)
plot(grid$strat)
plot(rasterToPolygons(grid$strat, dissolve = TRUE), col = "grey")

table(values(grid$depth))
table(values(grid$strat))

