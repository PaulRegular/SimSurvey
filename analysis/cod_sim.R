
library(SimSurvey)

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
                    export_dir = "analysis/cod_sim_exports/2018-10-26_age_clust_test")
# sim <- resume_test(export_dir = "analysis/cod_sim_exports/2018-10-26_age_clust_test")
setMKLthreads() # turn hyperthreading on again


# ## visualize results
# load("analysis/cod_sim_exports/2018-10-26_age_clust_test/test_output.RData")
# vis_sim(sim)


## Simulate same distribution across ages --------------------------------------

## Same as above cod-like simulation, except spatial distribution is the same
## across all ages
rm(pop)
rm(sim)
gc()

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
                                       length_group = 3,
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
                    export = "analysis/cod_sim_exports/2018-10-28_no_age_clust_test")
# sim <- resume_test(export_dir = "analysis/cod_sim_exports/2018-10-28_no_age_clust_test")
setMKLthreads() # turn hyperthreading on again


# ## visualize results
# load("analysis/cod_sim_exports/2018-10-28_no_age_clust_test/test_output.RData")
# vis_sim(sim)



## Test survey with strat-specific age-length-keys --------------------------------------

library(SimSurvey)

set.seed(438)
pop <- sim_abundance() %>%
  sim_distribution()


default_survey <- pop %>%
  sim_survey(n_sims = 10)

alt_survey <- pop %>%
  sim_survey(n_sims = 10, ages_cap = 5, age_length_group = 5,
             age_space_group = "strat")











