
library(SimSurvey)

## Simulate cod-like population
## See "imitate_cod_data.R" file for details on parameter choices

set.seed(438)
pop <- sim_abundance(ages = 1:20,
                     years = 1:20,
                     R = sim_R(mean = 30000000,
                               log_sd = 0.5,
                               random_walk = TRUE),
                     Z = sim_Z(mean = 0.5,
                               log_sd = 0.2,
                               phi_age = 0.9,
                               phi_year = 0.5),
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1,
                                       length_group = 3,
                                       digits = 0)) %>%
  sim_distribution(grid = sim_grid(x_range = c(-140, 140),
                                   y_range = c(-140, 140),
                                   res = c(3.5, 3.5),
                                   shelf_depth = 200,
                                   shelf_width = 100,
                                   depth_range = c(0, 1000),
                                   n_div = 1,
                                   strat_breaks = seq(0, 1000, by = 20),
                                   strat_splits = 2),
                   ays_covar = sim_ays_covar(sd = 2.8,
                                             range = 300,
                                             phi_age = 0.5,
                                             phi_year = 0.9,
                                             group_ages = 5:20),
                   depth_par = sim_parabola(mu = 200,
                                            sigma = 70))


## Test a series of surveys
## Simulate surveys and compare stratified estimates to the true index
setMKLthreads(1) # turn off MKL hyperthreading
surveys <- expand_surveys(set_den = c(0.5, 1, 2, 5, 10) / 1000,
                          lengths_cap = c(5, 10, 20, 50, 100, 500, 1000),
                          ages_cap = c(2, 5, 10, 20, 50))
surveys[surveys$set_den == 0.002 &
          surveys$lengths_cap == 500 &
          surveys$ages_cap == 10, ]    ## survey 98 ~ roughly current protocol
sim <- test_surveys(pop,
                    surveys = surveys,
                    keep_details = 98,
                    n_sims = 5,
                    n_loops = 200,
                    cores = 7,
                    q = sim_logistic(k = 2, x0 = 3),
                    export_dir = "analysis/cod_sim_exports/2018-09-17_age_clust_test")
# sim <- resume_test(export_dir = "analysis/cod_sim_exports/2018-09-17_age_clust_test")
setMKLthreads() # turn hyperthreading on again


# ## visualize results
# load("analysis/cod_sim_exports/2018-07-13_test/test_output.RData")
# vis_sim(sim)


# ## Exports for Geoff -----------------------------------------------------------
#
# library(data.table)
#
# ## ~ current sampling protocol
# base_case <- sim_survey_parallel(pop, n_sims = 5, n_loops = 5, cores = 6,
#                                  set_den = 2e-03, lengths_cap = 500, ages_cap = 10,
#                                  q = sim_logistic(k = 2, x0 = 3),
#                                  growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0),
#                                  quiet = FALSE)
# setdet <- base_case$setdet[, list(year, x, y, depth, strat, n, sim, set)]
# samp <- base_case$samp[base_case$samp$measured, ]
# samp$age[!samp$aged] <- NA
# samp$measured <- samp$aged <- NULL
# samp <- merge(setdet[, list(year, sim, set)], samp, by = "set")
# fwrite(setdet, file = "analysis/cod_sim_exports/2018-09-12_for_Geoff/01_setdet.csv")
# fwrite(samp, file = "analysis/cod_sim_exports/2018-09-12_for_Geoff/01_samp.csv")
# rm(base_case)
# gc()
#
# ## One example of a protocol that results in bias
# bias_case <- sim_survey_parallel(pop, n_sims = 1, n_loops = 25, cores = 6,
#                                  set_den = 1e-02, lengths_cap = 10, ages_cap = 10,
#                                  q = sim_logistic(k = 2, x0 = 3),
#                                  growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0),
#                                  quiet = FALSE)
# setdet <- bias_case$setdet[, list(year, x, y, depth, strat, n, sim, set)]
# samp <- bias_case$samp[bias_case$samp$measured, ]
# samp$age[!samp$aged] <- NA
# samp$measured <- samp$aged <- NULL
# samp <- merge(setdet[, list(year, sim, set)], samp, by = "set")
# fwrite(setdet, file = "analysis/cod_sim_exports/2018-09-12_for_Geoff/02_setdet.csv")
# fwrite(samp, file = "analysis/cod_sim_exports/2018-09-12_for_Geoff/02_samp.csv")
# rm(bias_case)
# gc()


## Simulate same distribution across ages --------------------------------------

## Same as above cod-like simulation, except spatial distribution is the same
## across all ages
rm(pop)
rm(sim)
gc()

set.seed(438)
pop <- sim_abundance(ages = 1:20,
                     years = 1:20,
                     R = sim_R(mean = 30000000,
                               log_sd = 0.5,
                               random_walk = TRUE),
                     Z = sim_Z(mean = 0.5,
                               log_sd = 0.2,
                               phi_age = 0.9,
                               phi_year = 0.5),
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1,
                                       length_group = 3,
                                       digits = 0)) %>%
  sim_distribution(grid = sim_grid(x_range = c(-140, 140),
                                   y_range = c(-140, 140),
                                   res = c(3.5, 3.5),
                                   shelf_depth = 200,
                                   shelf_width = 100,
                                   depth_range = c(0, 1000),
                                   n_div = 1,
                                   strat_breaks = seq(0, 1000, by = 20),
                                   strat_splits = 2),
                   ays_covar = sim_ays_covar(sd = 2.8,
                                             range = 300,
                                             phi_age = 0.5,
                                             phi_year = 0.9,
                                             group_ages = 1:20),
                   depth_par = sim_parabola(mu = 200,
                                            sigma = 70))

## Test a series of surveys
## Simulate surveys and compare stratified estimates to the true index
setMKLthreads(1) # turn off MKL hyperthreading
surveys <- expand_surveys(set_den = c(0.5, 1, 2, 5, 10) / 1000,
                          lengths_cap = c(5, 10, 20, 50, 100, 500, 1000),
                          ages_cap = c(2, 5, 10, 20, 50))
surveys[surveys$set_den == 0.002 &
          surveys$lengths_cap == 500 &
          surveys$ages_cap == 10, ]    ## survey 98 ~ roughly current protocol
sim <- test_surveys(pop,
                    surveys = surveys,
                    keep_details = 98,
                    n_sims = 5,
                    n_loops = 200,
                    cores = 7,
                    q = sim_logistic(k = 2, x0 = 3),
                    export = "analysis/cod_sim_exports/2018-09-18_no_age_clust_test")
# sim <- resume_test(dir = "analysis/cod_sim_exports/2018-09-18_no_age_clust_test")
setMKLthreads() # turn hyperthreading on again


# ## visualize results
# load("analysis/cod_sim_exports/2018-07-16_test/test_output.RData")
# vis_sim(sim)


