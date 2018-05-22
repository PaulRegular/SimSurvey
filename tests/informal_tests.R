
library(SimSurvey)

pop <- sim_abundance(years = 1:100)
vis_sim(pop)

pop <- sim_distribution(pop = sim_abundance(years = 1:30),
                        grid = sim_grid(res = c(3.5, 3.5)),
                        space_covar = sim_sp_covar(range = 50, sd = 0.1),
                        ay_covar = sim_ay_covar(sd = 10,
                                                phi_age = 0.5,
                                                phi_year = 0.5),
                        depth_par = sim_parabola(alpha = 0, sigma = 50))
vis_sim(pop)



pop <- sim_distribution(pop = sim_abundance(years = 1:5),
                        grid = sim_grid(res = c(3.5, 3.5)),
                        space_covar = sim_sp_covar(range = 50, sd = 0.1),
                        ay_covar = sim_ay_covar(sd = 10,
                                                phi_age = 0.5,
                                                phi_year = 0.5),
                        depth_par = sim_parabola(alpha = 0, sigma = 50))
vis_sim(pop)

sets <- sim_sets(pop)

