
library(SimSurvey)

sim <- sim_abundance(years = 1:100)
# vis_sim(sim)

set.seed(438)
sim <- sim_abundance(years = 1:20, R = sim_R(mean = 1e+09)) %>%
  sim_distribution(grid = sim_grid(res = c(3.5, 3.5)),
                   space_covar = sim_sp_covar(range = 50, sd = 0.1),
                   ay_covar = sim_ay_covar(sd = 10, phi_age = 0.5, phi_year = 0.5),
                   depth_par = sim_parabola(alpha = 0, sigma = 50)) %>%
  sim_survey(n_sims = 20)
# vis_sim(sim)

strat <- run_strat(sim)

library(plotly)
plot_ly() %>%
  add_lines(data = strat$strat2, x = ~year, y = ~total, split = ~sim,
            size = I(1), color = I("grey"), alpha = 0.4,
            showlegend = FALSE) %>%
  add_lines(x = sim$year, y = colSums(sim$I), color = I("black"), size = I(1)) %>%
  layout(yaxis = list(title = "I"),
         xaxis = list(title = "Year"))




