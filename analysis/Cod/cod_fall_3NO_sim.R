# Simulate Fall 3NO cod population
## See "imitate_cod_fall_3NO.R" file for details on parameter choices

library(SimSurvey)

set.seed(438)
pop <- sim_abundance(ages = 1:20,
                     years = 1:24,
                     R = sim_R(log_mean = log(75000000),
                               log_sd = 0.5,
                               random_walk = TRUE),
                     Z = sim_Z(log_mean = log(0.63),
                               log_sd = 0.3,
                               phi_age = 0.9,
                               phi_year = 0.5),
                     N0 = sim_N0(N0 = "exp", plot = FALSE),
                     growth = sim_vonB(Linf = 120, L0 = 5,
                                       K = 0.11, log_sd = 0.15,
                                       length_group = 3, digits = 0)) %>%
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
                   ays_covar = sim_ays_covar(sd = 2.5,
                                             range = 1700,
                                             lambda = .5,
                                             model = "matern",
                                             phi_age = 0.5,
                                             phi_year = 0.9,
                                             group_ages = 12:20),
                   depth_par = sim_parabola(mu = log(80),
                                            sigma = 0.25,
                                            sigma_right = 0.44, log_space = TRUE))


## Test a series of surveys
## Simulate surveys and compare stratified estimates to the true index
## Include baseline, increase and decrease of 80% and 50%; except for age, only 50%
surveys <- expand_surveys(set_den = c(0.5, 0.8, 1, 1.2, 1.5) / 1000,
                          lengths_cap = c(400, 500, 600),
                          ages_cap = c(5, 10, 15))
surveys[surveys$set_den == 0.001 &
          surveys$lengths_cap == 500 &
          surveys$ages_cap == 10, ]    ## survey 38 ~current protocol
sim <- test_surveys(pop,
                    surveys = surveys,
                    keep_details = 23,
                    n_sims = 5,
                    n_loops = 200,
                    cores = 7,
                    q = sim_logistic(k = 1.6, x0 = 1.7),
                    export_dir = "analysis/sim_exports/cod_sim_exports/2021-02-20_age_clust_test")

# sim <- resume_test(export_dir = "analysis/cod_sim_exports/2018-10-26_age_clust_test")


## Visualize Results
load("analysis/sim_exports/cod_sim_exports/2021-02-20_age_clust_test/test_output.RData")

plot_total_strat_fan(sim)
plot_length_strat_fan(sim, surveys = 1, years = 1:20, lengths = 0.5:144.5)
plot_age_strat_fan(sim, surveys = 1, years = 10, ages = 1:26)

plot_error_surface(sim)
plot_survey_rank(sim)

## Test survey with set-specific age-length-keys --------------------------------------

library(SimSurvey)

## Test alternate survey with strat specific age sampling and age-length-keys
set.seed(889)
surveys <- expand_surveys(set_den = c(1, 1.5) / 1000,
                          lengths_cap = c(250, 400, 500, 600, 750),
                          ages_cap = c(5, 10, 15))
sim <- test_surveys(pop,
                    surveys = surveys,
                    keep_details = 1,
                    n_sims = 5,
                    n_loops = 200,
                    cores = 7,
                    q = sim_logistic(k = 1.6, x0 = 1.7),
                    export_dir = "analysis/sim_exports/cod_sim_exports/2021-02-22_set_alk",
                    age_length_group = 2,
                    age_space_group = "set",
                    alk_scale = "set")


plot_age_strat_fan(sim, surveys = 5, ages = 1:10, years = 7)

plot_survey_rank(sim, which_strat = "length")
