

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
  sim_survey(n_sims = 10, age_sampling = "random") %>%
  run_strat() %>%
  strat_error()
rm(error)

## Test a series of surveys
## Simulate surveys and compare stratified estimates to the true index
setMKLthreads(1) # turn off MKL hyperthreading
surveys <- expand_surveys(set_den = c(2) / 1000,
                          lengths_cap = c(100),
                          ages_cap = c(3, 4))
sim <- test_surveys(pop,
                    surveys = surveys,
                    n_sims = 10,
                    n_loops = 100,
                    cores = 6,
                    q = sim_logistic(k = 2, x0 = 3),
                    age_sampling = "random") # export_dir = "tests/exports")
# sim <- resume_test(dir = "tests/exports")
setMKLthreads() # turn hyperthreading on again

plot_age_strat_fan(sim, surveys = 1:8, select_by = "age",
                   ages = sim$ages, years = sim$years)
sim$age_strat_error_stats
## Appears to be poorer performance with random vs length-stratified age sampling strategies



## Check consequences of narrow strata since it potentially contributes
## to the downward bias observed under the low set density scenario
library(raster)
grid <- sim_grid(x_range = c(-140, 140),
                 y_range = c(-140, 140),
                 res = c(3.5, 3.5),
                 shelf_depth = 250,
                 shelf_width = 0,
                 depth_range = c(0, 500),
                 n_div = 1,
                 strat_breaks = seq(0, 500, by = 25),
                 strat_splits = 2,
                 method = "linear")
plot(rasterToPolygons(grid$strat, dissolve = TRUE))
plot(grid)

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
  sim_distribution(grid = grid,
                   ays_covar = sim_ays_covar(sd = 2.8,
                                             range = 300,
                                             phi_age = 0.5,
                                             phi_year = 0.9,
                                             group_ages = 5:20),
                   depth_par = sim_parabola(mu = 200,
                                            sigma = 70))

error <- pop %>%
  sim_survey(n_sims = 10, age_sampling = "random") %>%
  run_strat() %>%
  strat_error()

plot_samp_dist(error, which_year = 4)

rm(error)


## Test a series of surveys
## Simulate surveys and compare stratified estimates to the true index
setMKLthreads(1) # turn off MKL hyperthreading
surveys <- expand_surveys(set_den = c(0.0005),
                          lengths_cap = c(100),
                          ages_cap = c(10))
sim <- test_surveys(pop,
                    surveys = surveys,
                    n_sims = 10,
                    n_loops = 100,
                    cores = 6,
                    q = sim_logistic(k = 2, x0 = 3)) # export_dir = "tests/exports")
# sim <- resume_test(dir = "tests/exports")
setMKLthreads() # turn hyperthreading on again

plot_total_strat_fan(sim)
plot_distribution_slider(sim, ages = 1:20, years = 6)
plot_samp_dist(sim, which_year = 6, which_sim = 5)
mean(sim$total_strat_error$error)
hist(sim$total_strat_error[year == 6]$error, breaks = 100, xlab = "error", main = "")
abline(v = 0, col = "red")
abline(v = mean(sim$total_strat_error$error), col = "blue")
abline(v = median(sim$total_strat_error$error), col = "green")
## Note that the percentiles in the fan plot make it look slightly biased...but the mean would look better
## because it lands in a higher place on a skewed distribution



