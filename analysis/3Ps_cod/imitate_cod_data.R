
## Explore cod sampling data from the RV survey and try and imitate
## the data using SimSurvey

## Survey grid -----------------------------------------------------------------

library(raster)
library(SimSurvey)

plot(survey_grid) # see data-raw folder
prod(res(survey_grid)) * sum(!is.na(values(survey_grid$cell)))
# 69147.59
length(unique(survey_grid$strat))
# 44
mean(table(values(survey_grid$strat)) * prod(res(survey_grid)))
# 1567.925
range(values(survey_grid$depth), na.rm = TRUE)
# -7.372526 920.745896

grid <- make_grid(x_range = c(-140, 140), y_range = c(-140, 140), res = c(3.5, 3.5),
                 shelf_depth = 200, shelf_width = 100, depth_range = c(0, 1000),
                 n_div = 1, strat_breaks = seq(0, 1000, by = 40), strat_splits = 2,
                 method = "spline")
prod(res(grid)) * ncell(grid)
# 78400
length(unique(grid$strat))
# 44
mean(table(values(grid$strat)) * prod(res(grid)))
# 1781.818
range(values(grid$depth), na.rm = TRUE)
# 15 940
plot(grid)
plot(grid$depth)
plot(grid$strat)
plot(rasterToPolygons(grid$strat, dissolve = TRUE), col = "grey")

## default settings of sim_grid() are similar to 3Ps


## Abundance and distribution --------------------------------------------------

library(Rstrap)
library(plotly)
library(data.table)

## REAL DATA

## Load survey data compiled for Rstrap
##load("analysis/rv_data/converted_set_details_2019-10-03.Rdata")
##load("analysis/rv_data/converted_length-frequency_data_2019-10-03.Rdata")
##load("analysis/rv_data/age-growth_data_2019-10-03.Rdata")

## Subset data to cod
con.setdet <- con.setdet[(con.setdet$rec == 5 | (con.setdet$rec == 6 & con.setdet$spec == 438)), ]
con.lf <- con.lf[con.lf$spec == 438, ]
ag <- ag[ag$spec == 438, ]
rv_data <- list(setdet = con.setdet, lf = con.lf, ag = ag)

## Save Rstrap output
index.strata <- c(293:300, 306:326, 705:708, 711:716, 779:783)
out <- strat.fun(setdet = rv_data$setdet, lf = rv_data$lf, ag = rv_data$ag,
                 data.series = "Campelen", program = "strat2 & strat1",
                 species = 438, survey.year = 1995:2017, season = "spring",
                 NAFOdiv = "3P", strat = c(293:300, 306:326, 705:708, 711:716, 779:783),
                 sex = c("male","female","unsexed"), length.group = 3,
                 group.by = "length & age",
                 export = NULL, plot.results = FALSE)

## Convert lat and lon to UTM
setdet <- data.table(out$raw.data$set.details)
xy <- cbind(-setdet$long.start, setdet$lat.start) %>%
  rgdal::project(., proj = proj4string(survey_grid)) %>%
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
set.seed(438)
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


## COMPS

## Compare set catch
hist(setdet$number, breaks = 100, xlab = "set catch", main = "real data")
hist(survey$setdet$n, breaks = 100, xlab = "set catch", main = "simulated data")

## Compare annual index
data_I <- out$strat2$abundance$summary$total[survey$years]
names(data_I) <- survey$years
sim_I <- colSums(survey$I)
barplot(data_I, names.arg = names(data_I), xlab = "age", main = "real data")
barplot(sim_I, names.arg = names(sim_I), xlab = "age", main = "simulated data")

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


