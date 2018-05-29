
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
vis_sim(abundance)


## Distribution ----------------------------------------------------------------

library(Rstrap)

## Load survey data compiled for Rstrap
load("analysis/rv_data/converted_set_details_2018-02-26.Rdata")
load("analysis/rv_data/converted_length-frequency_data_2018-02-26.Rdata")
load("analysis/rv_data/age-growth_data_2018-02-26.Rdata")

## Subset data to cod
con.setdet <- con.setdet[(con.setdet$rec == 5 | (con.setdet$rec == 6 & con.setdet$spec == 438)), ]
con.lf <- con.lf[con.lf$spec == 438, ]
setdet <- setdet[(setdet$rec == 5 | (setdet$rec == 6 & setdet$spec == 438)), ]
lf <- lf[lf$spec == 438, ]
ag <- ag[ag$spec == 438, ]
rv_data <- list(setdet = con.setdet, lf = con.lf, ag = ag)

## Save Rstrap output
index <- strat.fun(setdet = rv_data$setdet,
                   lf = rv_data$lf,
                   ag = rv_data$ag,
                   program = "strat2 & strat1",
                   data.series = c("Engel","Campelen"),
                   species = 438,
                   survey.year = 1983:2017,
                   season = "fall",
                   NAFOdiv = c("2J", "3K", "3L"),
                   strat = c(201:211, 213:217, 222, 223, 227:229, 234, 235, 237, 238, 240,
                             617:631, 633:640, 645, 650, 328, 341:350, 363:366, 368:372,
                             384:392),
                   sex = c("male", "female", "unsexed"),
                   length.group = 3,
                   group.by = "length & age",
                   plot.results = FALSE,
                   export = NULL)


setdet <- index$raw.data$set.details

