
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
mean(table(values(survey_grid$strat)) * res(survey_grid))
# 449.4992
range(values(survey_grid$depth), na.rm = TRUE)
# -7.372526 920.745896

grid <- sim_grid(x_range = c(-140, 140), y_range = c(-140, 140), res = c(3.5, 3.5),
                 shelf_depth = 200, shelf_width = 100, depth_range = c(0, 1000),
                 n_div = 1, strat_breaks = seq(0, 1000, by = 20), strat_splits = 2)
prod(res(grid)) * ncell(grid)
# 78400
length(unique(grid$strat))
# 42
mean(table(values(grid$strat)) * res(grid))
# 533.3333
range(values(grid$depth), na.rm = TRUE)
# 15 940
plot(grid)
plot(grid$depth)
plot(grid$strat)
plot(rasterToPolygons(grid$strat, dissolve = TRUE), col = "grey")

## default settings of sim_grid() are similar to 3Ps, except here we have 2 divisions
## could use 3Ps survey_grid, but will use a sim_grid for simplicity


## Abundance -------------------------------------------------------------------

## Roughly based on parameter estimates from NCAM
abundance <- sim_abundance(ages = 1:20, years = 1:20,
                           R = sim_R(mean = 1000000, log_sd = 0.4,
                                     random_walk = TRUE),
                           Z = sim_Z(mean = 0.7, log_sd = 0.2,
                                     phi_age = 0.9, phi_year = 0.5))
# vis_sim(abundance)


## Distribution ----------------------------------------------------------------

library(Rstrap)
library(plotly)

## Load survey data compiled for Rstrap
load("analysis/rv_data/converted_set_details_2018-02-26.Rdata")
load("analysis/rv_data/converted_length-frequency_data_2018-02-26.Rdata")
load("analysis/rv_data/age-growth_data_2018-02-26.Rdata")

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

## Melt age frequency data
af <- data.table::melt(setdet,
                       id.vars = c("survey.year", "vessel", "trip", "set", "easting", "northing"),
                       measure.vars = names(setdet)[grepl("^af", names(setdet))],
                       variable.name = "age", value.name = "freq")
af <- af[af$age != "afNA", ]
af$age <- as.integer(gsub("af", "", af$age))

## Hold age, animate year
plot_ly(data = af[af$age == 5, ]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~survey.year,
              sizes = c(5, 500), showlegend = FALSE) %>%
  animation_opts(frame = 5)

## Hold year, animate age
plot_ly(data = af[af$survey.year == 2009, ]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~age,
              sizes = c(5, 500), showlegend = FALSE) %>%
  animation_opts(frame = 5)

## Younger ages (< 4) are somewhat random (because of distribution or catchability?),
## hoever, correlation is strong across age. Less strong through time.

## Simulate data for comparison
abundance <- sim_abundance(ages = 1:20, years = 1:20,
                           R = sim_R(mean = 100000000, log_sd = 0.5,
                                     random_walk = TRUE),
                           Z = sim_Z(mean = 0.5, log_sd = 0.2,
                                     phi_age = 0.9, phi_year = 0.5))
grid <- sim_grid(x_range = c(-140, 140), y_range = c(-140, 140), res = c(3.5, 3.5),
                 shelf_depth = 200, shelf_width = 100, depth_range = c(0, 1000),
                 n_div = 1, strat_breaks = seq(0, 1000, by = 20), strat_splits = 2)
distribution <- sim_distribution(abundance, grid = grid,
                                 space_covar = sim_sp_covar(range = 30, sd = 0.3),
                                 ay_covar = sim_ay_covar(sd = 0.1,
                                                         phi_age = 0, phi_year = 0.2),
                                 depth_par = sim_parabola(mu = 250, sigma = 50))
survey <- sim_survey(distribution, n_sims = 1, light = FALSE,
                     set_den = 4 / 1000, lengths_cap = 400, ages_cap = 10,
                     q = sim_logistic(k = 2, x0 = 2.5),
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0))

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
points(length ~ age, pch = ".", data = survey$samp)
## Fish are caught based on age, not length...that's why there is a distinction
## between the two










sim_af <- survey$full_setdet

## Hold age, animate year
plot_ly(data = sim_af[sim_af$age == 5, ]) %>%
  add_markers(x = ~x, y = ~y, size = ~n, frame = ~year,
              sizes = c(5, 500), showlegend = FALSE) %>%
  animation_opts(frame = 5)

## Hold year, animate age
plot_ly(data = sim_af[sim_af$year == 10, ]) %>%
  add_markers(x = ~x, y = ~y, size = ~n, frame = ~age,
              sizes = c(5, 500), showlegend = FALSE) %>%
  animation_opts(frame = 5)

## Revisit the covar sd's and keep tweaking until you land on something that looks right



space_covar <- sim_sp_covar(range = 30, sd = 1)
xy <- survey$grid_xy[, c("x", "y")]
Sigma_space <- space_covar(xy)
w <- t(chol(Sigma_space))
nc <- nrow(xy)
xy$e <- w %*% rnorm(nc, 0, 1)
plot_ly(data = xy, x = ~x, y = ~y, z = ~e) %>% add_heatmap()






