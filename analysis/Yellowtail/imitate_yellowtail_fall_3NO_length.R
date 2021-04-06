# Yellowtail 3NO Fall
### Length based simulation ###

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
rv_data <- list(setdet = con.setdet, lf = con.lf)


## Investigate Rstrap output combined sexes
out <- strat.fun(setdet = rv_data$setdet, lf = rv_data$lf, ag = NULL,
                 data.series = "Campelen", program = "strat2 & strat1", which.survey = "multispecies",
                 species = 891, survey.year = c(1996:2019),
                 season = "fall",  # no age-growth 2002-2019, no data for 1995
                 NAFOdiv = c("3N", "3O"), strat = NULL,
                 sex = c("female", "male", "unsexed"),
                 length.group = 2, length.weight = NULL,
                 group.by = "length",
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

mean(setdet$set.depth.mean)


### Melt length frequency data
lf <- data.table::melt(setdet,
                       id.vars = c("survey.year", "vessel", "trip", "set", "easting", "northing", "set.depth.mean"),
                       measure.vars = names(setdet)[grepl("^lf", names(setdet))],
                       variable.name = "length", value.name = "freq")
lf <- lf[lf$length != "lfNA", ]
lf$length <- as.integer(gsub("lf", "", lf$length))


## Real by length
real_l <- data.frame(lf[lf$freq>0])
real_l  %>%
  ggplot(aes(x=set.depth.mean, y=freq, col=length))+
  geom_point() + scale_color_gradientn(colours = rainbow(5)) + theme_bw()

real_l  %>%
  ggplot(aes(x=set.depth.mean, y=freq)) +
  geom_point() + facet_wrap(~length)


## SIMULATED DATA

## Simulate data for comparison

# Abundance parameters and catchability curve roughly based on NCAM estimates
# Distribution parameters manually tweaked until results roughly corresponded
# to observations from 3NO plaice
set.seed(891)
pop <- sim_abundance(ages = 1:10,
                     years = 1:23,
                     R = sim_R(log_mean = log(7000000000),
                               log_sd = 0.6,
                               random_walk = FALSE),
                     Z = sim_Z(log_mean = log(0.64),
                               log_sd = 0.1,
                               phi_age = 0.2,
                               phi_year = 0.4),
                     N0 = sim_N0(N0 = "exp", plot = FALSE),
                     growth = sim_vonB(Linf = 56, L0 = 0,
                                       K = 0.13, log_sd = 0.12,
                                       length_group = 2, digits = 0)) %>%
  sim_distribution(grid,
                   ays_covar = sim_ays_covar(sd = 1.4,
                                             range = 190,
                                             phi_age = 0.8,
                                             phi_year = 0.7),
                                             #group_ages = 5:9),
                   depth_par = sim_parabola(mu = log(65),
                                            sigma = 0.15,
                                            sigma_right = 0.14, log_space = TRUE))


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

# Min sets per strata for all surveys are 2, set density is 2 for Yellowtail in 3NO,
# size of length bins for stratified age sampling is 2,
# max number of lengths per set are 300 and age group is 25
# All other values are default values (except q, which is modified based on annual index)

survey <- sim_survey(pop,
                     n_sims = 1,
                     q = sim_logistic(k = 1.6, x0 = 5.5),
                     trawl_dim = c(1.5, 0.02),
                     resample_cells = FALSE,
                     binom_error = TRUE,
                     min_sets = 2,
                     set_den = 1/1000,
                     lengths_cap = 300,
                     ages_cap = 25,
                     age_sampling = "stratified",
                     age_length_group = 2,
                     age_space_group = "division",
                     light = FALSE)

## COMPARISONS BETWEEN REAL AND SIMUALTED DATA

## Compare set catch

# Modified through sim_distribution arguments.

data_Z <- setdet[setdet$number==0,]
data_I <- setdet[setdet$number>0,]
sim_Z <- survey$setdet[survey$setdet$n==0,]
sim_I <- survey$setdet[survey$setdet$n>0,]
hist(data_Z$number, breaks = 100, xlab = "set catch", main = "set catch - real data")
hist(sim_Z$n, breaks = 100, xlab = "set catch", main = "set catch - simulated data")
hist(data_I$number, breaks = 100, xlab = "set catch", main = "set catch - real data")
hist(sim_I$n, breaks = 100, xlab = "set catch", main = "set catch - simulated data")

mean(data_I$number)
mean(sim_I$n)

plot_ly() %>%
  add_histogram(x = data_I$number, name = "real") %>%
  add_histogram(x = sim_I$n, name = "simulated") %>%
  layout(title = "Set catch")

## Compare annual index

data_I <- out$strat2$abundance$summary$total[survey$years]
names(data_I) <- survey$years
sim_I <- colSums(survey$I)
barplot(data_I, names.arg = names(data_I), xlab = "year", main = "annual index - real data")
barplot(sim_I, names.arg = names(sim_I), xlab = "year", main = "annual index - simulated data")

median(data_I)
mean(sim_I)

plot_ly() %>%
  add_lines(data = out$strat2$abundance$summary,
            x = ~seq_along(survey.year), y = ~total, name = "real") %>%
  add_lines(x = seq(survey$years), y = colSums(survey$I), name = "simulated") %>%
  layout(title = "Annual index", xaxis = list(title = "Year"))

## Compare index at length

data_I <- out$strat1$length$abundance$details
sim_I <- survey$samp[aged == TRUE, ]
nrow(data_I)
nrow(sim_I)
hist(data_I$length, xlab = "length", main = "real data", breaks = 200)
hist(sim_I$length, xlab = "length", main = "simulated data", breaks = 200)

mean(data_I$length)
mean(sim_I$length)

plot_ly(alpha = 0.9, nbinsx = 100) %>%
  add_histogram(x = data_I$length, name = "real") %>%
  add_histogram(x = sim_I$length, name = "simulated") %>%
  layout(title = "Average Index at Length", xaxis = list(title = "Length"))

##

data_I <- out$strat1$length$abundance$annual.totals
data_I_l <- rowMeans(data_I[data_I$length %in% survey$lengths, grepl("y", names(data_I))])
barplot(data_I_l, names.arg = data_I$length , xlab = "length", main = "Index @ length - real data")

sim_I_l <- rowMeans(survey$I_at_length)
barplot(sim_I_l, names.arg = names(sim_I_l), xlab = "length", main = " Index @ length- simulated data")


mean(data_I$length)
mean(sim_I$length)

plot_ly() %>%
  add_lines(x = seq_along(data_I_l), y = data_I_l, name = "real") %>%
  add_lines(x = seq_along(sim_I_l), y = sim_I_l, name = "simulated") %>%
  layout(title = "Average index at length", xaxis = list(title = "Length"))


## Compare age growth data
# Modify with sim_vonB (K) argument

data_I <- out$strat1$length$abundance$details
sim_I <- survey$samp[aged == TRUE, ]
nrow(data_I)
nrow(sim_I)
hist(data_I$length, xlab = "length", main = "age growth data - real data", breaks = 100)
hist(sim_I$length, xlab = "length", main = "age growth data - simulated data", breaks = 100)

mean(data_I$length)
mean(sim_I$length)

plot_ly(alpha = 0.6, nbinsx = 30) %>%
  add_histogram(x = data_I$length, name = "real") %>%
  add_histogram(x = sim_I$length, name = "simulated") %>%
  layout(title = "Average Index at Length", xaxis = list(title = "Length"))


## Compare relationship between catch and depth

data_I <- setdet[setdet$number>0,]
sim_I <- survey$setdet[survey$setdet$n>0,]
plot(as.numeric(data_I$set.depth.mean), data_I$number, xlab = "depth",
     ylab = "number", main = "real data", xlim = c(0, 1000))
plot(sim_I$depth, sim_I$n, xlab = "depth",
     ylab = "number", main = "simulated data", xlim = c(0, 1000))

median(data_I$set.depth.mean)
median(sim_I$depth)

plot_ly() %>%
  add_markers(x = data_I$set.depth.mean, y = data_I$number, name = "real") %>%
  add_markers(x = sim_I$depth, sim_I$n, name = "simulated") %>%
  layout(title = "Compare Catch Depth", xaxis = list(title = "Depth"),
         yaxis = list(title = "Number"))


## Compare intra-haul correlation
idvar <- c("vessel", "trip", "set", "year")
sub_lf <- merge(out$raw.data$set.details[, idvar], con.lf, by = idvar)
sub_lf$set <- as.numeric(as.factor(with(sub_lf, paste(vessel, trip, set, year))))
ind <- grep("^bin|^set$", names(con.lf)) # get length frequencies and expand to samples
lf <- sub_lf[, ind]
lf <- as.data.table(lf)
len_samp <- data.table::melt(lf, id.vars = "set")
len_samp <- len_samp[len_samp$value > 0, ]
len_samp$value <- round(len_samp$value)
len_samp$length <- as.numeric(gsub("bin", "", len_samp$variable))
len_samp <- as.data.frame(len_samp)
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
        inches = 0.1, main = "size of distribution - real data",
        xlab = "x", ylab = "y")
symbols(survey$setdet$x, survey$setdet$y,
        circles = sqrt(survey$setdet$n / pi),
        inches = 0.1, main = "size of distribution - simulated data",
        xlab = "x", ylab = "y")



## Real data (hold length or year and animate the other)
plot_ly(data = lf[lf$length == 30, ]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~survey.year,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 5)
plot_ly(data = lf[lf$survey.year == 2000,]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~age,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)

# Examine at the age dimension, with frequency scaled by age to allow for
# distribution shifts at older ages to be visible
lf[lf$survey.year == 1998, ] %>%
  group_by(length) %>%
  mutate(scaled_freq = scale(freq)) %>%
  plot_ly() %>%
  add_markers(x = ~easting, y = ~northing, size = ~scaled_freq, frame = ~length,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)

## Real data all lengths and years
plot_ly(data = lf) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~survey.year,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 5)
plot_ly(data = lf) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~length,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)

## Again, scale within age
lf %>%
  group_by(length) %>%
  mutate(scaled_freq = scale(freq)) %>%
  plot_ly() %>%
  add_markers(x = ~easting, y = ~northing, size = ~scaled_freq, frame = ~length,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)

## Simulated data

sim_lf <- data.frame(survey$full_setdet)
sim_lf %>%
  filter(length == 20) %>%
  group_by(year) %>%
  plot_ly(x = ~x, y = ~y, size = ~n, frame = ~year,
          sizes = c(5, 1000), showlegend = FALSE) %>%
  add_markers() %>%
  animation_opts(frame = 5)

sim_lf %>%
  filter(year == 5) %>%
  group_by(length) %>%
  plot_ly(x = ~x, y = ~y, size = ~n, frame = ~age,
          sizes = c(5, 1000), showlegend = FALSE) %>%
  add_markers() %>%
  animation_opts(frame = 500)


