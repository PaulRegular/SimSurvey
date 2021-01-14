# Simulate Fall 3NO female plaice population
## See "imitate_plaice_fall_3NO_female.R" file for details on parameter choices

library(SimSurvey)

set.seed(889)
pop <- sim_abundance(ages = 1:19,
                     years = 1:20,
                     R = sim_R(log_mean = log(410000000),
                               log_sd = 0.9,
                               random_walk = FALSE),
                     Z = sim_Z(log_mean = log(0.44),
                               log_sd = 0.6,
                               phi_age = 0.5,
                               phi_year = 0.5),
                     N0 = sim_N0(N0 = "exp", plot = FALSE),
                     growth = sim_vonB(Linf = 50.58, L0 = 3,  # Fitted for male growth
                                       K = 0.12, log_sd = 0.1,
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
                   ays_covar = sim_ays_covar(sd = 2.7,
                                             range = 900,
                                             phi_age = 0.9,
                                             phi_year = 0.5,
                                             group_ages = 16:19),
                   depth_par = sim_parabola(mu = log(75),
                                            sigma = 0.15,
                                            sigma_right = 0.5, log_space = TRUE))


## Test a series of surveys
## Simulate surveys and compare stratified estimates to the true index
## Include baseline, increase and decrease of 80% and 50%; except for age, only 50%
surveys <- expand_surveys(set_den = c(0.5, 0.8, 1, 1.2, 1.5) / 1000,
                          lengths_cap = c(75, 120, 150, 180, 225),
                          ages_cap = c(10, 20, 30))
surveys[surveys$set_den == 0.001 &
          surveys$lengths_cap == 150 &
          surveys$ages_cap == 20, ]    ## survey 38 ~current protocol
sim <- test_surveys(pop,
                    surveys = surveys,
                    keep_details = 38,
                    n_sims = 5,
                    n_loops = 200,
                    cores = 7,
                    q = sim_logistic(k = 1.4, x0 = 3.3),
                    export_dir = "analysis/plaice_male_sim_exports/2021-01-11_age_clust_test")
# sim <- resume_test(export_dir = "analysis/cod_sim_exports/2018-10-26_age_clust_test")


## Visualize Results
load("analysis/plaice_male_sim_exports/2021-01-11_age_clust_test/test_output.RData")

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
                          lengths_cap = c(75, 120, 150, 180, 225),
                          ages_cap = c(10, 20, 30))
sim <- test_surveys(pop,
                    surveys = surveys,
                    keep_details = 1,
                    n_sims = 5,
                    n_loops = 200,
                    cores = 7,
                    q = sim_logistic(k = 1.4, x0 = 3.3),
                    export_dir = "analysis/plaice_male_sim_exports/2021-01-11_set_alk",
                    age_length_group = 2,
                    age_space_group = "set",
                    alk_scale = "set")


plot_age_strat_fan(sim, surveys = 5, ages = 1:10, years = 7)

plot_survey_rank(sim, which_strat = "length")
