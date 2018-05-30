
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

grid <- sim_grid(x_range = c(-150, 150), y_range = c(-150, 150), res = c(3.5, 3.5),
              depth_range = c(1, 500), n_div = 2, strat_breaks = seq(0, 500, by = 20))
prod(res(grid)) * ncell(grid)
# 90601
length(unique(grid$strat))
# 50
mean(table(values(grid$strat)) * res(grid))
# 517.72
plot(grid)

## default settings of sim_grid() are similar to 3Ps, except here we have 2 divisions
## could use 3Ps survey_grid, but will use a sim_grid for simplicity


## Abundance -------------------------------------------------------------------

## Roughly based on parameter estimates from NCAM
abundance <- sim_abundance(years = 1:10,
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
                 export = NULL)

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

## Roughly based on parameter estimates from NCAM
grid <- sim_grid(x_range = c(-150, 150), y_range = c(-150, 150), res = c(5, 5),
                 depth_range = c(1, 500), n_div = 2, strat_breaks = seq(0, 500, by = 20))
abundance <- sim_abundance(years = 1:10,
                           R = sim_R(mean = 100000000, log_sd = 0.4,
                                     random_walk = TRUE),
                           Z = sim_Z(mean = 0.5, log_sd = 0.2,
                                     phi_age = 0.9, phi_year = 0.5))
distribution <- sim_distribution(abundance, grid = grid,
                                 space_covar = sim_sp_covar(range = 50, sd = 0.1),
                                 ay_covar = sim_ay_covar(sd = 10, phi_age = 0.8, phi_year = 0),
                                 depth_par = sim_parabola(alpha = 0, mu = 250, sigma = 50))
survey <- sim_survey(distribution, n_sims = 1, light = FALSE)

sim_af <- survey$full_setdet

## Hold age, animate year
plot_ly(data = sim_af[sim_af$age == 5, ]) %>%
  add_markers(x = ~x, y = ~y, size = ~n, frame = ~year,
              sizes = c(5, 500), showlegend = FALSE) %>%
  animation_opts(frame = 5)

## Hold year, animate age
plot_ly(data = sim_af[sim_af$year == 1, ]) %>%
  add_markers(x = ~x, y = ~y, size = ~n, frame = ~age,
              sizes = c(5, 500), showlegend = FALSE) %>%
  animation_opts(frame = 5)
