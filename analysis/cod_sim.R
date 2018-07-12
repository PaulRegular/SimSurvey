
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
                   ays_covar = sim_ays_covar(sd = 2.8,
                                             range = 300,
                                             phi_age = 0.1,
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
res <- test_surveys(pop,
                    surveys = surveys,
                    keep_details = 98,
                    n_sims = 5,
                    n_loops = 200,
                    cores = 6,
                    q = sim_logistic(k = 2, x0 = 3),
                    growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0),
                    export = "analysis/cod_sim_exports")
# res <- resume_test(dir = "analysis/cod_sim_exports/2018-07-07_test")
setMKLthreads() # turn hyperthreading on again


## visualize results
load("analysis/cod_sim_exports/2018-07-07_test/test_output.RData")
vis_sim(sim)

abun_plot <- plot_total_abundance(sim)
set_plot <- plot_set_catch(sim)
export(abun_plot, file = "analysis/paper/figures/total_abundance.png",
       vwidth = 1000, vheight = 800, zoom = 20)

# plotly_IMAGE(abun_plot, out_file = "analysis/paper/figures/total_abundance.png")




