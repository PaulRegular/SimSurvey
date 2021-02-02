## Make grid function for 3L

## Survey grid -----------------------------------------------------------------

library(sf)
library(raster)
library(SimSurvey)
library(dplyr)
library(ggplot2)
library(plotly)

## UTM projection for Newfoundland
utm_proj <- "+proj=utm +zone=21 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

## Import strata shapefile
strat_polys <- st_read("data-raw/DFO_NL_survey_strat/2HJ3KLNOP_Strata_Polygons_WGS84.shp",
                       layer = "2HJ3KLNOP_Strata_Polygons_WGS84")

## Subset to Division
strat_polys_div <- strat_polys %>%
  filter (DIV == "3L")

strat_polys_div_utm <- strat_polys_div %>% st_transform(utm_proj)

## Import bathy data
strat_bathy <- raster("data-raw/GEBCO_derived_bathy/2HJ3KLNOP_GEBCO_bathy.nc")
strat_bathy <- raster::mask(strat_bathy, strat_polys_div) # cells within strata

## Strata-by-strata depths
strat_no <- unique(strat_polys_div$STRAT)
means <- numeric(length(strat_no))
mins <- numeric(length(strat_no))
maxs <- numeric(length(strat_no))
names(means) <- names(mins) <- names(maxs) <- strat_no

for (s in strat_no) {

  temp <- raster::mask(strat_bathy, strat_polys_div[strat_polys_div$STRAT == s, ])

  cat(s, "\n")
  cat("mean:", mean(temp[], na.rm = TRUE), "\n")
  cat("min:", min(temp[], na.rm = TRUE), "\n")
  cat("max:", max(temp[], na.rm = TRUE), "\n")

  means[as.character(s)] <- mean(temp[], na.rm = TRUE)
  mins[as.character(s)] <- min(temp[], na.rm = TRUE)
  maxs[as.character(s)] <- max(temp[], na.rm = TRUE)

}

depth_by_strata <- data.frame(strat = strat_no,
                              mean = means,
                              min = mins, max = maxs)

## Survey area
survey_area <- sum(sf::st_area(strat_polys_div_utm))
survey_area
# 165955.1 [km^2]

## Area of the strata
strat_sums <- strat_polys_div_utm %>%
  mutate(area = sf::st_area(strat_polys_div_utm)) %>%
  group_by(STRAT) %>%
  summarize(strat_area = sum(area))
range(strat_sums$strat_area)
# Units: [km^2]
# 219.8599 9815.2453

## Number of strata
strata <- length(unique(strat_polys_div_utm$STRAT))
strata
# 69

## Mean area
mean_area <- mean(strat_sums$strat_area)
mean_area
# 2405.146 [km^2]

## Range depth in the survey area
depth_range <- range(strat_bathy[], na.rm = TRUE)
depth_range
# -1851    19

## Mean depth in the survey area
depth_mean <- mean(strat_bathy[], na.rm = TRUE)
depth_mean
# Units: [km^2]
# -261.4167

depth_median <- median(strat_bathy[], na.rm = TRUE)
depth_median
# Units: [km^2]
# -161

## Determines x_y_range
sqrt(survey_area)/2
# 203.6879 [km]

## TO DO: Modify make_grid arguments to create a similar survey area
## shelf depth/width modifies mean grid depth, start breaks/splits modifies
## strata (length), method determines slope of shelf

grid <- make_grid(x_range = c(-204, 204),
                  y_range = c(-204, 204),
                  res = c(3.5, 3.5),
                  shelf_depth = 200,
                  shelf_width = 100,
                  depth_range = c(0, 1850),
                  n_div = 1,
                  strat_breaks = seq(0, 2000, by = 50),
                  strat_splits = 2,
                  method = "bezier")

tibble::lst(survey_area, strata, mean_area, depth_range, depth_mean, depth_median)

prod(res(grid)) * ncell(grid)
#
length(unique(grid$strat))
#
mean(table(values(grid$strat)) * prod(res(grid)))
#
range(values(grid$depth), na.rm = TRUE)
#
mean(values(grid$depth), na.rm = TRUE)
#
median(values(grid$depth), na.rm = TRUE)
#

xyz <- data.frame(rasterToPoints(grid$depth))
plot_ly(data = xyz, x = ~x, y = ~-depth) %>% add_lines()

plot(grid)
plot(grid$depth)
plot(grid$strat)
plot(rasterToPolygons(grid$strat, dissolve = TRUE), col = "grey")

grid_dat <- data.frame(raster::rasterToPoints(grid))

strat_depths <- grid_dat %>%
  group_by(strat) %>%
  summarize(mean = mean(depth), min = min(depth), max = max(depth))

## Compare real vs simulated depths

real_depths <- data.frame(depth_by_strata)
real_depths <- abs(real_depths)
sim_depths <- data.frame(strat_depths)
all_depths <- rbind(real_depths, sim_depths)

all_depths <- all_depths %>% mutate(grid = factor(ifelse(strat > 300, "real", "sim")))

all_depths %>%
  filter(!is.na(grid)) %>%
  ggplot(aes(x=mean, color=grid, fill=grid)) +
  geom_histogram(position="dodge", bins=10, alpha = 1) + theme_bw()

plot_grid(grid)
