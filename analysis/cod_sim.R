
library(SimSurvey)

## Simulate cod-like population
## See "imitate_cod_data.R" file for details on parameter choices

set.seed(438)
pop <- sim_abundance(ages = 1:20,
                     years = 1:20,
                     R = sim_R(mean = 100000000,
                               log_sd = 0.5,
                               random_walk = TRUE),
                     Z = sim_Z(mean = 0.5,
                               log_sd = 0.2,
                               phi_age = 0.9,
                               phi_year = 0.5)) %>%
  sim_distribution(grid = sim_grid(x_range = c(-140, 140),
                                   y_range = c(-140, 140),
                                   res = c(3.5, 3.5),
                                   shelf_depth = 200,
                                   shelf_width = 100,
                                   depth_range = c(0, 1000),
                                   n_div = 1,
                                   strat_breaks = seq(0, 1000, by = 20),
                                   strat_splits = 2),
                   space_covar = sim_sp_covar(range = 40,
                                              sd = 0.1),
                   ay_covar = sim_ay_covar(sd = 10,
                                           phi_age = 0.1,
                                           phi_year = 0.8,
                                           group_ages = 5:20),
                   depth_par = sim_parabola(mu = 250,
                                            sigma = 50))

## Test a series of surveys
## Simulate surveys and compare stratified estimates to the true index
setMKLthreads(1) # turn off MKL hyperthreading
res <- test_surveys(pop,
                    # surveys = expand_surveys(set_den = c(0.3, 0.5, 0.8,
                    #                                      1, 2, 3, 6, 9) / 1000,
                    #                          lengths_cap = c(2, 3, 5, 8, 10, 20,
                    #                                          30, 60, 90, 100, 200,
                    #                                          400, 600, 1000),
                    #                          ages_cap = c(2, 3, 5, 8, 10,
                    #                                       20, 30, 60)),
                    surveys = expand_surveys(set_den = c(0.5, 1, 3, 9) / 1000,
                                             lengths_cap = c(3, 8, 20, 60, 100,
                                                             400, 1000),
                                             ages_cap = c(3, 8, 20, 60)),
                    n_sims = 5,
                    n_loops = 50,
                    cores = 6,
                    q = sim_logistic(k = 2, x0 = 3),
                    growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0),
                    export = "analysis/cod_sim_exports")
# res <- resume_test(dir = "analysis/2018-05-31_test")
setMKLthreads() # turn hyperthreading on again




## Simulate one survey
survey <- sim_survey(pop,
                     n_sims = 1,
                     light = FALSE,
                     set_den = 3 / 1000,
                     lengths_cap = 400,
                     ages_cap = 10,
                     q = sim_logistic(k = 2, x0 = 3),
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0))





