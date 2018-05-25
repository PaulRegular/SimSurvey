
library(SimSurvey)

sim <- sim_abundance(years = 1:100)
# vis_sim(sim)

set.seed(438)
sim <- sim_abundance(years = 1:10, R = sim_R(mean = 1e+09)) %>%
  sim_distribution(grid = sim_grid(res = c(3.5, 3.5)),
                   space_covar = sim_sp_covar(range = 50, sd = 0.1),
                   ay_covar = sim_ay_covar(sd = 10, phi_age = 0.5, phi_year = 0.5),
                   depth_par = sim_parabola(alpha = 0, sigma = 50)) %>%
  sim_survey(n_sims = 1)
# vis_sim(sim)

data <- strat_data(sim)
strat <- run_strat(sim)

I <- sim$setdet[, list(total = sum(I)), by = "year"]$total
I_hat <- strat$strat2$total
plot(I_hat)
lines(I)
I_hat
I


