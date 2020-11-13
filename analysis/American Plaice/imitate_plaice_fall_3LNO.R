
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
  geom_histogram(position="dodge", bins=10, alpha = 1)+ theme_bw()

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
st_xy <- st_as_sf(data.frame(long = -setdet$long.start, lat = setdet$lat.start),
                  coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs(strat_polys_utm)) %>%
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

## ~ 0.005 sets per km^2
## ~ 3 sets per 200 sq. NM

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
                     Z = sim_Z(log_mean = log(0.5),   # natural (M:0.30) and fishing (F:0.20) mortality
                               log_sd = 0.2,
                               phi_age = 0.9,         # M decreases with age, F increases with age
                               phi_year = 0.5),       # inverted u-shape graph ages 9-14 between 1996-2017)
                     N0 = sim_N0(N0 = "exp", plot = FALSE),
                     growth = sim_vonB(Linf = 50, L0 = 3,       # Fitted for female growth
                                       K = 0.18, log_sd = 0.1,
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


## COMPS

## Compare set catch
hist(setdet$number, breaks = 100, xlab = "set catch", main = "real data")
hist(survey$setdet$n, breaks = 100, xlab = "set catch", main = "simulated data")

## Compare annual index
data_I <- out$strat2$abundance$summary$total[survey$years]
names(data_I) <- survey$years
sim_I <- colSums(survey$I)
barplot(data_I, names.arg = names(data_I), xlab = "age", main = "annual index - real data")
barplot(sim_I, names.arg = names(sim_I), xlab = "age", main = "annual - index - simulated data")

## Compare index at age
data_I <- out$strat1$age$abundance$annual.totals
data_I <- rowMeans(data_I[data_I$age %in% survey$ages, grepl("y", names(data_I))])
sim_I <- rowMeans(survey$I)
barplot(data_I, names.arg = names(data_I), xlab = "age", main = "real data")
barplot(sim_I, names.arg = names(sim_I), xlab = "age", main = "simulated data")

## Compare age growth data
data_I <- out$raw.data$age.growth
sim_I <- survey$samp[aged == TRUE, ]
nrow(data_I)
nrow(sim_I)
hist(data_I$length, xlab = "length", main = "real data", breaks = 100)
hist(sim_I$length, xlab = "length", main = "simulated data", breaks = 100)
plot(length ~ age, data = data_I, main = "real data",
     xlim = range(survey$ages))
plot(length ~ age, data = sim_I, main = "simulated data",
     xlim = range(survey$ages))
# points(length ~ age, pch = ".", data = survey$samp)
## Fish are caught based on age, not length...that's why there is a distinction
## between the two. If the catchability simulation were length based, a scattered
## older individual would be in the mix along the tail of the length distribution

## Compare relationship between catch and depth
## (could use to improve the accuracy of depth in the real data, but the pattern is clear)
data_I <- setdet
sim_I <- survey$setdet
plot(as.numeric(data_I$set.depth.min), data_I$number, xlab = "depth",
     ylab = "number", main = "real data", xlim = c(0, 1000))
plot(sim_I$depth, sim_I$n, xlab = "depth",
     ylab = "number", main = "simulated data", xlim = c(0, 1000))

## Compare intra-haul correlation
idvar <- c("vessel", "trip", "set", "year")
sub_lf <- merge(out$raw.data$set.details[, idvar], con.lf, by = idvar)
sub_lf$set <- as.numeric(as.factor(with(sub_lf, paste(vessel, trip, set, year))))
ind <- grep("^bin|^set$", names(con.lf)) # get length frequencies and expand to samples
lf <- sub_lf[, ind]
len_samp <- data.table::melt(lf, id.vars = "set")
len_samp <- len_samp[len_samp$value > 0, ]
len_samp$value <- round(len_samp$value)
len_samp$length <- as.numeric(gsub("bin", "", len_samp$variable))
len_samp <- len_samp[rep(row.names(len_samp), len_samp$value), c("set", "length")]
len_samp <- data.table(len_samp)
sub_sets <- sample(len_samp$set, size = 10)
stripchart(length ~ set, data = len_samp[set %in% sub_sets, ],
           vertical = TRUE, pch = 1, cex = 0.5, method = "jitter", jitter = 0.2,
           main = "real data")
icc(len_samp$length, len_samp$set)

len_samp <- survey$samp[survey$samp$measured, ]
sub_sets <- sample(len_samp$set, size = 10)
stripchart(length ~ set, data = len_samp[set %in% sub_sets, ],
           vertical = TRUE, pch = 1, cex = 0.5, method = "jitter", jitter = 0.2,
           main = "simulated data")
icc(len_samp$length, len_samp$set)



## Now size up the distribution

symbols(setdet$long.start, setdet$lat.start,
        circles = sqrt(setdet$number / pi),
        inches = 0.1, main = "real data",
        xlab = "x", ylab = "y")
symbols(survey$setdet$x, survey$setdet$y,
        circles = sqrt(survey$setdet$n / pi),
        inches = 0.1, main = "simulated data",
        xlab = "x", ylab = "y")


## Real data (hold age or year and animate the other)
plot_ly(data = af[af$age == 5, ]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~survey.year,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 5)
plot_ly(data = af[af$survey.year == 2013, ]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~age,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)
## Younger ages (< 4) are somewhat random (because of distribution or catchability?),
## hoever, correlation is strong across age. Less strong through time.


## Hold age or year and animate the other
# sim_af <- data.frame(survey$full_setdet)
# for (a in rev(survey$ages)) {
#   d <- sim_af[sim_af$year == 1 & sim_af$age == a, ]
#   radius <- sqrt( d$n / pi )
#   symbols(d$x, d$y, circles = radius, inches = 0.1, main = paste("age", a),
#           xlab = "x", ylab = "y")
# }
# for (y in rev(survey$years)) {
#   d <- sim_af[sim_af$year == y & sim_af$age == 1, ]
#   radius <- sqrt( d$n / pi )
#   symbols(d$x, d$y, circles = radius, inches = 0.1, main = paste("year", y),
#           xlab = "x", ylab = "y")
# }

sim_af <- data.frame(survey$full_setdet)
sim_af %>%
  filter(age == 5) %>%
  group_by(year) %>%
  plot_ly(x = ~x, y = ~y, size = ~n, frame = ~year,
          sizes = c(5, 1000), showlegend = FALSE) %>%
  add_markers() %>%
  animation_opts(frame = 5)
sim_af %>%
  filter(year == 5) %>%
  group_by(age) %>%
  plot_ly(x = ~x, y = ~y, size = ~n, frame = ~age,
          sizes = c(5, 1000), showlegend = FALSE) %>%
  add_markers() %>%
  animation_opts(frame = 500)
