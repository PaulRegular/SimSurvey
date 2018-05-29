
library(SimSurvey)

set.seed(438)
pop <- sim_abundance(years = 1:10, R = sim_R(mean = 1e+09)) %>%
  sim_distribution(grid = sim_grid(res = c(3.5, 3.5)),
                   space_covar = sim_sp_covar(range = 50, sd = 0.1),
                   ay_covar = sim_ay_covar(sd = 10, phi_age = 0.5, phi_year = 0.5),
                   depth_par = sim_parabola(alpha = 0, sigma = 50))
# vis_sim(pop)

survey <- sim_survey(pop, n_sims = 10) %>%
  run_strat() %>% strat_error()

setMKLthreads(1) # turn off MKL hyperthreading
res <- test_surveys(pop, n_sims = 5, n_loops = 100, cores = 8)
setMKLthreads() # turn hyperthreading on again


library(plotly)
plot_ly() %>%
  add_lines(data = res$total_strat_error, x = ~year, y = ~I_hat, split = ~sim,
            size = I(1), color = ~I("grey"), alpha = 0.4, frame = ~survey,
            showlegend = FALSE) %>%
  add_lines(x = ~unique(year), y = ~unique(I), color = I("black"), size = I(1)) %>%
  layout(yaxis = list(title = "I"),
         xaxis = list(title = "Year"))


