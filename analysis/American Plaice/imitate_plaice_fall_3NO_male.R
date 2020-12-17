# Plaice 3NO Fall Male

## Make grid function for each Division

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

# SUBSET TO DIVISION

strat_polys_div <- strat_polys %>%
                   filter (DIV == "3N" | DIV == "3O")

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

#x_y_range
sqrt(survey_area)/2
# 183.842 [km]

## TO DO: Modify make_grid arguments to create a similar survey area

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


## Abundance and distribution --------------------------------------------------

library(Rstrap)
library(plotly)
library(data.table)
library(tidyr)

## REAL DATA

## Subset data to plaice
con.setdet <- con.setdet[(con.setdet$rec == 5 | (con.setdet$rec == 6 & con.setdet$spec == 889)), ]
con.lf <- con.lf[con.lf$spec == 889, ]
ag <- ag[ag$spec == 889, ]
rv_data <- list(setdet = con.setdet, lf = con.lf, ag = ag)

## Investigate Rstrap output for both males and females as plaice
## have sexually dimorphic growth
out_all <- strat.fun(setdet = rv_data$setdet, lf = rv_data$lf, ag = rv_data$ag,
                 data.series = "Campelen", program = "strat2 & strat1", which.survey = "multispecies",
                 species = 889, survey.year = c(1995:2013, 2015:2017),
                 season = "fall",  # no age-growth 2014, 2018, 2019
                 NAFOdiv = c("3N", "3O"), strat = NULL,
                 sex = c("female", "male", "unsexed"), sep.sex = TRUE, # sexually dimorphic growth
                 indi.alk = c(TRUE,TRUE,FALSE),
                 length.group = 2, length.weight = NULL,
                 group.by = "length & age",
                 export = NULL, plot.results = FALSE)

## Examine sex ratio by AGE
out_all$strat1$age$abundance$summary %>%
  filter(sex %in% c("male", "female")) %>%
  group_by(age, sex) %>%
  summarise(total = sum(total)) %>%
  ungroup() %>%
  plot_ly(x = ~age, y = ~total, color = ~sex) %>%
  add_lines()

out_all$strat1$age$abundance$summary %>%
  filter(age == 16, sex %in% c("male", "female")) %>%
  plot_ly(x = ~survey.year, y = ~total, color = ~sex) %>%
  add_lines()

## Examine sex ratio by YEAR
out_all$strat1$age$abundance$summary %>%
  filter(sex %in% c("male", "female")) %>%
  group_by(survey.year, sex) %>%
  summarise(total = sum(total)) %>%
  ungroup() %>%
  plot_ly(x = ~survey.year, y = ~total, color = ~sex) %>%
  add_lines()

out_all$strat1$age$abundance$summary %>%
  filter(survey.year == 2016, sex %in% c("male", "female")) %>%
  plot_ly(x = ~age, y = ~total, color = ~sex) %>%
  add_lines()

## Examines sex ratio above or below 50%

out_all$strat1$age$abundance$summary %>%
  filter(sex %in% c("male", "female")) %>%
  select(survey.year, age, sex, total) %>%
  pivot_wider(names_from = sex, values_from = total) %>%
  mutate(ratio = male / (male + female)) %>%
  group_by(age) %>%
  summarise(mean_ratio = mean(ratio, na.rm = TRUE))

out_all$strat1$age$abundance$summary %>%
  filter(sex %in% c("male", "female")) %>%
  select(survey.year, age, sex, total) %>%
  pivot_wider(names_from = sex, values_from = total) %>%
  mutate(ratio = male / (male + female)) %>%
  filter(age < 10) %>%
  plot_ly(x = ~survey.year, y = ~ratio - 0.5, color = ~factor(age)) %>%
  add_lines()

## Males are greater than 50% abundance in the population before age 7,
## females are greater than 50% at age 8 and older
## Higher male recruitment and mortality than females needed in simulated survey

## Save Rstrap output by extracting only male data
out <- strat.fun(setdet = rv_data$setdet, lf = rv_data$lf, ag = rv_data$ag,
                     data.series = "Campelen", program = "strat2 & strat1", which.survey = "multispecies",
                     species = 889, survey.year = c(1995:2013, 2015:2017),
                     season = "fall",  # no age-growth 2014, 2018, 2019
                     NAFOdiv = c("3N", "3O"), strat = NULL,
                     sex = "male", sep.sex = FALSE,
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


## SIMULATED DATA

## Simulate data for comparison
## - Abundance parameters and catchability curve roughly based on NCAM estimates
## - Distribution parameters manually tweaked until results roughly corresponded
##   to observations from 3NO plaice
set.seed(889)
pop <- sim_abundance(ages = 1:19,
                     years = 1:20,
                     R = sim_R(log_mean = log(410000000),
                               log_sd = 0.9,
                               random_walk = FALSE),
                     Z = sim_Z(log_mean = log(0.44),
                               log_sd = 0.6,
                               phi_age = 0.5,
                               phi_year = 0.5),
                     N0 = sim_N0(N0 = "exp", plot = FALSE),
                     growth = sim_vonB(Linf = 50.58, L0 = 3,  # Fitted for male growth
                                       K = 0.12, log_sd = 0.1,
                                       length_group = 2, digits = 0)) %>%
  sim_distribution(grid,
                   ays_covar = sim_ays_covar(sd = 2.7,
                                             range = 900,
                                             phi_age = 0.9,
                                             phi_year = 0.5,
                                             group_ages = 16:19),
                   depth_par = sim_parabola(mu = log(75),
                                            sigma = 0.15,
                                            sigma_right = 0.5, log_space = TRUE))


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
                     q = sim_logistic(k = 1.4, x0 = 3.3),
                     trawl_dim = c(1.5, 0.02),
                     resample_cells = FALSE,
                     binom_error = TRUE,
                     min_sets = 2,
                     set_den = 1/1000,
                     lengths_cap = 150,
                     ages_cap = 20,
                     age_sampling = "stratified",
                     age_length_group = 2,
                     age_space_group = "division",
                     light = FALSE)


## COMPS

## Compare set catch
## FREQUENCY OF SIM DATA SHOULD BE HALF
data_Z <- setdet[setdet$number==0,]
data_I <- setdet[setdet$number>0,]
sim_Z <- survey$setdet[survey$setdet$n==0,]
sim_I <- survey$setdet[survey$setdet$n>0,]
hist(data_Z$number, breaks = 100, xlab = "set catch", main = "set catch - real data")
hist(sim_Z$n, breaks = 100, xlab = "set catch", main = "set catch - simulated data")
hist(data_I$number, breaks = 100, xlab = "set catch", main = "set catch - real data")
hist(sim_I$n, breaks = 100, xlab = "set catch", main = "set catch - simulated data")

median(data_I$number * 0.5)
median(sim_I$n)

## Note adjustment for sex ratio - current simulation is focused on females
## while the set details 'number' from the survey includes both male and females
plot_ly() %>%
  add_histogram(x = data_I$number * 0.5, name = "real") %>%
  add_histogram(x = sim_I$n, name = "simulated") %>%
  layout(title = "Set catch")


## Compare annual index
## FREQUENCY OF SIM DATA SHOULD BE HALF
data_I <- out$strat2$abundance$summary$total[survey$years]
names(data_I) <- survey$years
sim_I <- colSums(survey$I)
barplot(data_I, names.arg = names(data_I), xlab = "year", main = "annual index - real data")
barplot(sim_I, names.arg = names(sim_I), xlab = "year", main = "annual index - simulated data")

mean(data_I*.05)
mean(sim_I)

## Again, note the 50% adjustment of the 'strat2' totals since they are based
## on set-level numbers that include male and females
plot_ly() %>%
  add_lines(data = out$strat2$abundance$summary,
            x = ~seq_along(survey.year), y = ~total * 0.5, name = "real") %>%
  add_lines(x = seq(survey$years), y = colSums(survey$I), name = "simulated") %>%
  layout(title = "Annual index", xaxis = list(title = "Year"))


## Compare index at age
data_I <- out$strat1$age$abundance$annual.totals
data_I <- rowMeans(data_I[data_I$age %in% survey$ages, grepl("y", names(data_I))])
sim_I <- rowMeans(survey$I)
barplot(data_I, names.arg = names(data_I), xlab = "age", main = "index at age - real data")
barplot(sim_I, names.arg = names(sim_I), xlab = "age", main = "index at age - simulated data")

mean(data_I)
mean(sim_I)

## These values, in contrast to the 'strat2' values, are limited to females
## given the strat.fun call used above
plot_ly() %>%
  add_lines(x = seq_along(data_I), y = data_I, name = "real") %>%
  add_lines(x = seq_along(sim_I), y = sim_I, name = "simulated") %>%
  layout(title = "Average index at age", xaxis = list(title = "Age"))


## Compare age growth data
data_I <- out$raw.data$age.growth
sim_I <- survey$samp[aged == TRUE, ]
nrow(data_I)
nrow(sim_I)
hist(data_I$length, xlab = "length", main = "age growth data - real data", breaks = 100)
hist(sim_I$length, xlab = "length", main = "age growth data - simulated data", breaks = 100)

mean(data_I$length)
mean(sim_I$length)

plot_ly() %>%
  add_histogram(x = data_I$length, name = "real") %>%
  add_histogram(x = sim_I$length, name = "simulated") %>%
  layout(title = "Age Growth")

plot(length ~ age, data = data_I, main = "age growth data - real data",
     xlim = range(survey$ages))
plot(length ~ age, data = sim_I, main = "age growth data - simulated data",
     xlim = range(survey$ages))
# points(length ~ age, pch = ".", data = survey$samp)
## Fish are caught based on age, not length...that's why there is a distinction
## between the two. If the catchability simulation were length based, a scattered
## older individual would be in the mix along the tail of the length distribution

plot_ly() %>%
  add_markers(x = data_I$age - 0.25, y = data_I$length, name = "real") %>%
  add_markers(x = sim_I$age + 0.25, y = sim_I$length, name = "simulated") %>%
  layout(title = "Age Growth",
         xaxis = list(title = "Age"),
         yaxis = list(title = "Length"))



## Compare relationship between catch and depth
## (could use to improve the accuracy of depth in the real data, but the pattern is clear)
## FREQUENCY OF SIM DATA SHOULD BE HALF
data_I <- setdet[setdet$number>0,]
sim_I <- survey$setdet[survey$setdet$n>0,]
plot(as.numeric(data_I$set.depth.mean), data_I$number, xlab = "depth",
     ylab = "number", main = "real data", xlim = c(0, 1000))
plot(sim_I$depth, sim_I$n, xlab = "depth",
     ylab = "number", main = "simulated data", xlim = c(0, 1000))

mean(data_I$set.depth.mean)
mean(sim_I$depth)

## Again, note 0.5 adjustment for the set-level number
plot_ly() %>%
  add_markers(x = data_I$set.depth.mean, y = data_I$number * 0.5, name = "real") %>%
  add_markers(x = sim_I$depth, sim_I$n, name = "simulated") %>%
  layout(title = "Compare Catch Depth", xaxis = list(title = "Depth"),
                                                     yaxis = list(title = "Number"))


## Relationship of catch and depth by age

## Real by age
real_a <- data.frame(af[af$freq>0])
real_a  %>%
  ggplot(aes(x=set.depth.mean, y=freq, col=age))+
  geom_point() + scale_color_gradientn(colours = rainbow(5)) + theme_bw()

real_a  %>%
  ggplot(aes(x=set.depth.mean, y=freq)) +
  geom_point() + facet_wrap(~age)

## Real by agegroup
real_a <- real_a %>% mutate(agegroup = case_when(age %in% 1:19 ~ "age 1-19",
                                                 age %in% 20:25 ~ "age 20-25"))

real_a$agegroup <- as.factor(real_a$agegroup)
real_a %>% filter(!is.na(agegroup)) %>%
            filter(agegroup == "age 1-19") %>%
          ggplot(aes(x=set.depth.mean, y=freq, col=agegroup))+
          geom_point() + scale_color_brewer(palette="Spectral")


## Simulated
sim_a <- data.frame(survey$full_setdet[survey$full_setdet$n>0])
sim_a %>% ggplot(aes(x=depth, y=n,col=age)) +
  geom_point() + scale_color_gradientn(colours = rainbow(5)) + theme_bw()

sim_a <- sim_a %>% mutate(agegroup = case_when(age %in% 1:19 ~ "age 1-19",
                                               age %in% 20:26 ~ "age 20-26"))
sim_a$agegroup <- as.factor(sim_a$agegroup)
sim_a %>% filter(!is.na(agegroup)) %>%
  filter(agegroup == "age 1-19") %>%
  ggplot(aes(x=depth, y=n,col=agegroup)) +
  geom_point() +scale_color_brewer(palette="Spectral") + theme_bw()

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


## Real data (hold age or year and animate the other)
plot_ly(data = af[af$age == 7, ]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~survey.year,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 5)
plot_ly(data = af[af$survey.year == 1999,]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~age,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)

## Have a look at the age dimension, with freq scaled by age to allow dist
## shifts at older ages to be visible
af[af$survey.year == 2011, ] %>%
  group_by(age) %>%
  mutate(scaled_freq = scale(freq)) %>%
  plot_ly() %>%
  add_markers(x = ~easting, y = ~northing, size = ~scaled_freq, frame = ~age,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)


## Real data all ages and years
plot_ly(data = af) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~survey.year,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 5)
plot_ly(data = af) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~age,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)


## Again, scale within age
af %>%
  group_by(age) %>%
  mutate(scaled_freq = scale(freq)) %>%
  plot_ly() %>%
  add_markers(x = ~easting, y = ~northing, size = ~scaled_freq, frame = ~age,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)



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
