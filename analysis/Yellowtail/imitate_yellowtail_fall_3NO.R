# Yellowtail 3NO Spring

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
  filter (DIV == "3N" | DIV == "3O")

strat_polys_div_utm <- strat_polys_div %>% st_transform(utm_proj)

## Import bathy data
strat_bathy <- raster("data-raw/GEBCO_derived_bathy/2HJ3KLNOP_GEBCO_bathy.nc")
strat_bathy <- raster::mask(strat_bathy, strat_polys_div) # cells within strata


## Survey area
survey_area <- sum(sf::st_area(strat_polys_div_utm))
survey_area
# 135191.5 [km^2]

## Area of the strata
strat_sums <- strat_polys_div_utm %>%
  mutate(area = sf::st_area(strat_polys_div_utm)) %>%
  group_by(STRAT) %>%
  summarize(strat_area = sum(area))
range(strat_sums$strat_area)
# Units: [km^2]
# 207.7644 10359.0056

## Number of strata
strata <- length(unique(strat_polys_div_utm$STRAT))
strata
# 71

## Mean area
mean_area <- mean(strat_sums$strat_area)
mean_area
# 1904.106 [km^2]

## Range depth in the survey area
depth_range <- range(strat_bathy[], na.rm = TRUE)
depth_range
# -2102    -2

## Mean depth in the survey area
depth_mean <- mean(strat_bathy[], na.rm = TRUE)
depth_mean
# Units: [km^2]
# -206.2931

depth_median <- median(strat_bathy[], na.rm = TRUE)
depth_median
# Units: [km^2]
# -78

## Determines x_y_range
sqrt(survey_area)/2
# 183.842 [km]

## TO DO: Modify make_grid arguments to create a similar survey area
## shelf depth/width modifies mean grid depth, start breaks/splits modifies
## strata (length), method determines slope of shelf

grid <- make_grid(x_range = c(-184, 184),
                  y_range = c(-184, 184),
                  res = c(3.5, 3.5),
                  shelf_depth = 60,
                  shelf_width = 170,
                  depth_range = c(0, 1600),
                  n_div = 2,
                  strat_breaks = seq(0, 1600, by = 65),
                  strat_splits = 4,
                  method = "bezier")

tibble::lst(survey_area, strata, mean_area, depth_range, depth_mean, depth_median)

prod(res(grid)) * ncell(grid)
# 135056.2
length(unique(grid$strat))
# 72
mean(table(values(grid$strat)) * prod(res(grid)))
# 1875.781
range(values(grid$depth), na.rm = TRUE)
# 6 1415
mean(values(grid$depth), na.rm = TRUE)
# 208.0762
median(values(grid$depth), na.rm = TRUE)
# 83

xyz <- data.frame(rasterToPoints(grid$depth))
plot_ly(data = xyz, x = ~x, y = ~-depth) %>% add_lines()

plot(grid)
plot(grid$depth)
plot(grid$strat)
plot(rasterToPolygons(grid$strat, dissolve = TRUE), col = "grey")


## Abundance and distribution --------------------------------------------------

library(Rstrap)
library(plotly)
library(data.table)
library(tidyr)

## REAL DATA

## Subset data to plaice
con.setdet <- con.setdet[(con.setdet$rec == 5 | (con.setdet$rec == 6 & con.setdet$spec == 891)), ]
con.lf <- con.lf[con.lf$spec == 891, ]
ag <- ag[ag$spec == 891, ]
rv_data <- list(setdet = con.setdet, lf = con.lf, ag = ag)

## Investigate Rstrap output for both males and females
out_all <- strat.fun(setdet = rv_data$setdet, lf = rv_data$lf, ag = rv_data$ag,
                     data.series = "Campelen", program = "strat2 & strat1", which.survey = "multispecies",
                     species = 891, survey.year = c(1995:2001),
                     season = "fall",  # no age-growth 2002-2019
                     NAFOdiv = c("3N", "3O"), strat = NULL,
                     sex = c("female", "male", "unsexed"), sep.sex = TRUE, # if sexually dimorphic growth
                     indi.alk = c(TRUE,TRUE,FALSE),
                     length.group = 2, length.weight = NULL,
                     group.by = "length & age",
                     export = NULL, plot.results = FALSE)

## Examines sex ratio by AGE
out_all$strat1$age$abundance$summary %>%
  filter(sex %in% c("male", "female")) %>%
  group_by(age, sex) %>%
  summarise(total = sum(total)) %>%
  ungroup() %>%
  plot_ly(x = ~age, y = ~total, color = ~sex) %>%
  add_lines()

out_all$strat1$age$abundance$summary %>%
  filter(age == 1, sex %in% c("male", "female")) %>%
  plot_ly(x = ~survey.year, y = ~total, color = ~sex) %>%
  add_lines()

## Examines sex ratio by YEAR
out_all$strat1$age$abundance$summary %>%
  filter(sex %in% c("male", "female")) %>%
  group_by(survey.year, sex) %>%
  summarise(total = sum(total)) %>%
  ungroup() %>%
  plot_ly(x = ~survey.year, y = ~total, color = ~sex) %>%
  add_lines()

out_all$strat1$age$abundance$summary %>%
  filter(survey.year == 2000, sex %in% c("male", "female")) %>%
  plot_ly(x = ~age, y = ~total, color = ~sex) %>%
  add_lines()

## Examines sex ratio above or below 50%
out_all$strat1$age$abundance$summary %>%
  filter(sex %in% c("female", "male")) %>%
  select(survey.year, age, sex, total) %>%
  pivot_wider(names_from = sex, values_from = total) %>%
  mutate(ratio = female / (female + male)) %>%
  group_by(age) %>%
  summarise(mean_ratio = mean(ratio, na.rm = TRUE)) %>%
  print (n=25)

out_all$strat1$age$abundance$summary %>%
  filter(sex %in% c("male", "female")) %>%
  select(survey.year, age, sex, total) %>%
  pivot_wider(names_from = sex, values_from = total) %>%
  mutate(ratio = male / (male + female)) %>%
  filter(age < 10) %>%
  plot_ly(x = ~survey.year, y = ~ratio - 0.5, color = ~factor(age)) %>%
  add_lines()


## Investigate Rstrap output combined sexes
out <- strat.fun(setdet = rv_data$setdet, lf = rv_data$lf, ag = rv_data$ag,
                 data.series = "Campelen", program = "strat2 & strat1", which.survey = "multispecies",
                 species = 891, survey.year = c(1996:2019),
                 season = "fall",  # no age-growth 2002-2019, no data for 1995
                 NAFOdiv = c("3N", "3O"), strat = NULL,
                 sex = c("female", "male", "unsexed"),
                 length.group = 2, length.weight = NULL,
                 group.by = "length & age",
                 export = NULL, plot.results = FALSE)

## Convert lat and lon to UTM
setdet <- data.table(out$raw.data$set.details)
st_xy <- st_as_sf(data.frame(long = -setdet$long.start, lat = setdet$lat.start),
                  coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs(strat_polys_div_utm)) %>%
  st_coordinates() %>% data.frame(.)
names(st_xy) <- c("easting", "northing")
setdet <- cbind(setdet, st_xy)

## Set density
den <- setdet[rec == 5, list(n = .N), by = c("survey.year", "NAFOdiv", "strat", "strat.area")]
den <- den[, list(strat_area_nm = sum(strat.area),
                  strat_area_km = sum(strat.area * 3.4299),
                  den_km = sum(n) / sum(strat.area * 3.4299),
                  den_nm = sum(n) / sum(strat.area / 200)), by = c("survey.year")]
den

## ~ 0.001 sets per km^2
## ~ 1 sets per 200 sq. NM

## Melt age frequency data
af <- data.table::melt(setdet,
                       id.vars = c("survey.year", "vessel", "trip", "set", "easting", "northing", "set.depth.mean"),
                       measure.vars = names(setdet)[grepl("^af", names(setdet))],
                       variable.name = "age", value.name = "freq")
af <- af[af$age != "afNA", ]
af$age <- as.integer(gsub("af", "", af$age))
