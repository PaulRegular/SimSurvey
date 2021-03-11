# Simulate Fall 3NO yellowtail population
## See "imitate_yellowtail_fall_3NO.R" file for details on parameter choices

library(SimSurvey)

set.seed(891)
pop <- sim_abundance(ages = 1:10,
                     years = 1:6,
                     R = sim_R(log_mean = log(12000000000),
                               log_sd = 0.6,
                               random_walk = TRUE),
                     Z = sim_Z(log_mean = log(0.64),
                               log_sd = 0.1,
                               phi_age = 0.2,
                               phi_year = 0.4),
                     N0 = sim_N0(N0 = "exp", plot = FALSE),
                     growth = sim_vonB(Linf = 56, L0 = 0,
                                       K = 0.13, log_sd = 0.12,
                                       length_group = 2, digits = 0)) %>%
  sim_distribution(grid = make_grid(x_range = c(-184, 184),
                                    y_range = c(-184, 184),
                                    res = c(3.5, 3.5),
                                    shelf_depth = 60,
                                    shelf_width = 170,
                                    depth_range = c(0, 1600),
                                    n_div = 2,
                                    strat_breaks = seq(0, 1600, by = 65),
                                    strat_splits = 4,
                                    method = "bezier"),
                   ays_covar = sim_ays_covar(sd = 1.4,
                                             range = 190,
                                             phi_age = 0.8,
                                             phi_year = 0.7),
                   #group_ages = 5:9),
                   depth_par = sim_parabola(mu = log(65),
                                            sigma = 0.15,
                                            sigma_right = 0.14, log_space = TRUE))

## Test a series of surveys
## Simulate surveys and compare stratified estimates to the true index
## Include baseline, increase and decrease of 80% and 50%; except for age, only 50%
surveys <- expand_surveys(set_den = c(0.5, 0.8, 1, 1.2, 1.5) / 1000,
                          lengths_cap = c(150, 240, 300, 360, 450),
                          ages_cap = c(15, 25, 35))
surveys[surveys$set_den == 0.001 &
          surveys$lengths_cap == 300 &
          surveys$ages_cap == 25, ]    ## survey 38 ~current protocol
sim <- test_surveys(pop,
                    surveys = surveys,
                    keep_details = 38,
                    n_sims = 5,
                    n_loops = 200,
                    cores = 7,
                    q = sim_logistic(k = 1.6, x0 = 5.5),
                    export_dir = "analysis/sim_exports/yellowtail_sim_exports/2021-03-_age_clust_test")


## Visualize Results
load("analysis/sim_exports/yellowtail_sim_exports/2021-03-_age_clust_test/test_output.RData")

plot_total_strat_fan(sim)
plot_length_strat_fan(sim, surveys = 1, years = 1:20, lengths = 0.5:144.5)
plot_age_strat_fan(sim, surveys = 1, years = 10, ages = 1:26)

plot_error_surface(sim)
plot_survey_rank(sim)

## Test survey with set-specific age-length-keys --------------------------------------

library(SimSurvey)

## Test alternate survey with strat specific age sampling and age-length-keys
set.seed(891)
surveys <- expand_surveys(set_den = c(1, 1.2) / 1000,
                          lengths_cap = c(150, 240, 300, 360, 450),
                          ages_cap = c(15, 25, 35))
sim <- test_surveys(pop,
                    surveys = surveys,
                    keep_details = 1,
                    n_sims = 5,
                    n_loops = 200,
                    cores = 7,
                    q = sim_logistic(k = 1.6, x0 = 5.5),
                    export_dir = "analysis/sim_exports/yellowtail_sim_exports/2021-03-_set_alk",
                    age_length_group = 2,
                    age_space_group = "set",
                    alk_scale = "set")


plot_age_strat_fan(sim, surveys = 5, ages = 1:10, years = 7)

plot_survey_rank(sim, which_strat = "length")

