
## Explore American plaice sampling data for the fall season in Div 3LNO
## from the RV survey and try to imitate the data using SimSurvey

## Survey grid -----------------------------------------------------------------

library(sf)
library(raster)
library(SimSurvey)
library(dplyr)
library(ggplot2)

## UTM projection for Newfoundland
utm_proj <- "+proj=utm +zone=21 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

## Import 3LNO strata shapefile
strat_polys <- st_read("data-raw/DFO_NL_survey_strat/2HJ3KLNOP_Strata_Polygons_WGS84.shp",
                       layer = "2HJ3KLNOP_Strata_Polygons_WGS84")
strat_polys <- strat_polys[strat_polys$DIV %in% c("3L", "3N", "3O"), ] # subset to 3LNO
strat_polys_utm <- st_transform(strat_polys, utm_proj)

## Import bathy data
strat_bathy <- raster("data-raw/GEBCO_derived_bathy/2HJ3KLNOP_GEBCO_bathy.nc")
strat_bathy <- raster::mask(strat_bathy, strat_polys) # cells within strata


## Strata-by-strata depths
strat_no <- unique(strat_polys$STRAT)
means <- numeric(length(strat_no))
mins <- numeric(length(strat_no))
maxs <- numeric(length(strat_no))
names(means) <- names(mins) <- names(maxs) <- strat_no

for (s in strat_no) {

  temp <- raster::mask(strat_bathy, strat_polys[strat_polys$STRAT == s, ])

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
survey_area <- sum(sf::st_area(strat_polys_utm))
survey_area
# 301146.6 [km^2]

## Area of the strata
strat_sums <- strat_polys_utm %>%
  mutate(area = sf::st_area(strat_polys_utm)) %>%
  group_by(STRAT) %>%
  summarize(strat_area = sum(area))
range(strat_sums$strat_area)
# Units: [km^2]
# 207.7644 10359.0056

## Number of strata
strata <- length(unique(strat_polys_utm$STRAT))
strata
# 140

## Mean area
mean_area <- mean(strat_sums$strat_area)
mean_area
# 2151.047 [km^2]

## Range depth in the survey area
depth_range <- range(strat_bathy[], na.rm = TRUE)
depth_range
# -2102    19

## Mean depth in the survey area
depth_mean <- mean(strat_bathy[], na.rm = TRUE)
depth_mean
# Units: [km^2]
# -237.2556

#x_y_range
sqrt(survey_area)/2
# 275

## TO DO: Modify make_grid arguments to create a similar survey area

grid <- make_grid(x_range = c(-275, 275),
                  y_range = c(-275, 275),
                  res = c(3.5, 3.5),
                  shelf_depth = 200,
                  shelf_width = 225,
                  depth_range = c(0, 2100),
                  n_div = 3,
                  strat_breaks = seq(0, 1000, by = 45),
                  strat_splits = 3,
                  method = "spline")

tibble::lst(survey_area, strata, mean_area, depth_range, depth_mean)

prod(res(grid)) * ncell(grid)
# 301950.2
length(unique(grid$strat))
# 140
mean(table(values(grid$strat)) * prod(res(grid)))
# 2074.363
range(values(grid$depth), na.rm = TRUE)
# 1 1972
mean(values(grid$depth), na.rm = TRUE)
# 244.2038

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
  geom_histogram(position="dodge", bins=10, alpha = 1)

plot_grid(grid)


## Abundance and distribution --------------------------------------------------

library(Rstrap)
library(plotly)
library(data.table)

## REAL DATA

## Subset data to plaice
con.setdet <- con.setdet[(con.setdet$rec == 5 | (con.setdet$rec == 6 & con.setdet$spec == 889)), ]
con.lf <- con.lf[con.lf$spec == 889, ]
ag <- ag[ag$spec == 889, ]
rv_data <- list(setdet = con.setdet, lf = con.lf, ag = ag)

## Save Rstrap output
out <- strat.fun(setdet = rv_data$setdet, lf = rv_data$lf, ag = rv_data$ag,
                 data.series = "Campelen", program = "strat2 & strat1", which.survey = "multispecies",
                 species = 889, survey.year = c(1995:2003, 2005:2013, 2015:2017),
                 season = "fall",  # incomplete coverage of 3L in 2004, 2014; no age-growth 2018, 2019
                 NAFOdiv = c("3L", "3N", "3O"), strat = NULL,
                 sex = c("male","female","unsexed"), sep.sex = TRUE, # sexually dimorphic growth
                 indi.alk = c(TRUE,TRUE,FALSE),
                 length.group = 2, length.weight = NULL,
                 group.by = "length & age",
                 export = NULL, plot.results = FALSE)

## Convert lat and lon to UTM
setdet <- data.table(out$raw.data$set.details)
xy <- cbind(-setdet$long.start, setdet$lat.start) %>%
  rgdal::project(., proj = crs(strat_polys_utm)) %>%
  data.frame(.)
names(xy) <- c("easting", "northing")
setdet <- cbind(setdet, xy)

## Set density
den <- setdet[rec == 5, list(n = .N), by = c("survey.year", "strat", "strat.area")]
den <- den[, list(strat_area_nm = sum(strat.area),
                  strat_area_km = sum(strat.area * 3.4299),
                  den_km = sum(n) / sum(strat.area * 3.4299),
                  den_nm = sum(n) / sum(strat.area / 200)), by = c("survey.year")]
den
## ~ 0.003 sets per km^2
## ~ 2 sets per 200 sq. NM

## Melt age frequency data
af <- data.table::melt(setdet,
                       id.vars = c("survey.year", "vessel", "trip", "set", "easting", "northing"),
                       measure.vars = names(setdet)[grepl("^af", names(setdet))],
                       variable.name = "age", value.name = "freq")
af <- af[af$age != "afNA", ]
af$age <- as.integer(gsub("af", "", af$age))


## SIMULATED DATA

## Simulate data for comparison
## - Abundance parameters and catchability curve roughly based on NCAM estimates
## - Distribution parameters manually tweaked until results roughly corresponded to
##   observations from 3Ps cod
set.seed(889)
pop <- sim_abundance(ages = 1:20,
                     years = 1:20,
                     R = sim_R(log_mean = log(30000000),
                               log_sd = 0.5,
                               random_walk = TRUE),
                     Z = sim_Z(log_mean = log(0.5),
                               log_sd = 0.2,
                               phi_age = 0.9,
                               phi_year = 0.5),
                     #N0 = sim_N0(N0 = "exp", plot = FALSE),
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, log_sd = 0.1,
                                       length_group = 3, digits = 0)) %>%
  sim_distribution(grid = make_grid(x_range = c(-140, 140),
                                    y_range = c(-140, 140),
                                    res = c(3.5, 3.5),
                                    shelf_depth = 200,
                                    shelf_width = 100,
                                    depth_range = c(0, 1000),
                                    n_div = 1,
                                    strat_breaks = seq(0, 1000, by = 40),
                                    strat_splits = 2,
                                    method = "spline"),
                   ays_covar = sim_ays_covar(sd = 2.8,
                                             range = 300,
                                             phi_age = 0.5,
                                             phi_year = 0.9,
                                             group_ages = 5:20),
                   depth_par = sim_parabola(mu = 200,
                                            sigma = 70))

## Quick look at distribution
sp_N <- data.frame(merge(pop$sp_N, pop$grid_xy, by = "cell"))
for (j in rev(pop$ages)) {
  z <- xtabs(N ~ x + y, subset = age == j & year == 1, data = sp_N)
  image(z = z, axes = FALSE, col = viridis::viridis(100), main = paste("age", j))
}
for (i in rev(pop$years)) {
  z <- xtabs(N ~ x + y, subset = age == 1 & year == i, data = sp_N)
  image(z = z, axes = FALSE, col = viridis::viridis(100), main = paste("year", i))
}


survey <- sim_survey(pop,
                     n_sims = 1,
                     light = FALSE,
                     set_den = 3 / 1000,
                     lengths_cap = 400,
                     ages_cap = 10,
                     q = sim_logistic(k = 2, x0 = 3))

