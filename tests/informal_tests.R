

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

error <- pop %>%
  sim_survey(n_sims = 10) %>%
  run_strat() %>%
  strat_error()
rm(error)

## Test a series of surveys
## Simulate surveys and compare stratified estimates to the true index
setMKLthreads(1) # turn off MKL hyperthreading
surveys <- expand_surveys(set_den = c(10) / 1000,
                          lengths_cap = c(10, 1000),
                          ages_cap = c(10))
sim <- test_surveys(pop,
                    surveys = surveys,
                    n_sims = 5,
                    n_loops = 100,
                    cores = 7,
                    q = sim_logistic(k = 2, x0 = 3))
setMKLthreads() # turn hyperthreading on again




