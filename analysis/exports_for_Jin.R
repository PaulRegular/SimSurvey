
library(SimSurvey)
library(data.table)

## Simulate cod-like population
## See "imitate_cod_data.R" file for details on parameter choices

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
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1,
                                       log_sd = 0.1, length_group = 3,
                                       digits = 0)) %>%
  sim_distribution(grid = make_grid(x_range = c(-140, 140),
                                    y_range = c(-140, 140),
                                    res = c(3.5, 3.5),
                                    shelf_depth = 200,
                                    shelf_width = 100,
                                    depth_range = c(0, 1000),
                                    n_div = 1,
                                    strat_breaks = seq(0, 1000, by = 40),
                                    strat_splits = 2),
                   ays_covar = sim_ays_covar(sd = 2.8,
                                             range = 300,
                                             phi_age = 0.5,
                                             phi_year = 0.9,
                                             group_ages = 5:20),
                   depth_par = sim_parabola(mu = 200,
                                            sigma = 70))



## ~ current sampling protocol
base_case <- sim_survey_parallel(pop, n_sims = 5, n_loops = 5, cores = 6,
                                 set_den = 2e-03, lengths_cap = 500, ages_cap = 10,
                                 q = sim_logistic(k = 2, x0 = 3),
                                 quiet = FALSE)

setdet <- base_case$setdet
samp <- base_case$samp
fwrite(setdet, file = "analysis/cod_sim_exports/2020-02-06_for_Jin/sim_cod_3Ps_setdet.csv")
fwrite(samp, file = "analysis/cod_sim_exports/2020-02-06_for_Jin/sim_cod_3Ps_samp.csv")
save(base_case, file = "analysis/cod_sim_exports/2020-02-06_for_Jin/sim_for_Jin.RData")

rm(base_case)
gc()



## 2020-02-06 - real 3Ps data --------------------------------------------------

library(Rstrap)

## Load survey data compiled for Rstrap
load("analysis/rv_data/converted_set_details_2019-10-03.Rdata")
load("analysis/rv_data/converted_length-frequency_data_2019-10-03.Rdata")
load("analysis/rv_data/age-growth_data_2019-10-03.Rdata")

## Subset data to cod
con.setdet <- con.setdet[(con.setdet$rec == 5 |
                            (con.setdet$rec == 6 & con.setdet$spec == 438)) &
                           con.setdet$NAFOdiv == "3P", ]
con.lf <- con.lf[con.lf$spec == 438 & con.lf$NAFOdiv == "3P", ]
ag <- ag[ag$spec == 438 & ag$NAFOdiv == "3P", ]
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

fwrite(out$raw.data$set.details,
       file = "analysis/cod_sim_exports/2020-02-06_for_Jin/real_cod_3Ps_setdet.csv")
fwrite(out$raw.data$age.growth,
       file = "analysis/cod_sim_exports/2020-02-06_for_Jin/real_cod_3Ps_samp.csv")


